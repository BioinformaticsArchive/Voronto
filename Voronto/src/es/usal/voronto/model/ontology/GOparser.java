package es.usal.voronto.model.ontology;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.emory.mathcs.backport.java.util.Arrays;
import es.usal.voronto.model.expression.ExpressionData;
import es.usal.voronto.view.VoronoiVisualization;

import keggapi.Definition;
import keggapi.KEGGLocator;
import keggapi.KEGGPortType;

public class GOparser {
	
	public static void main(String[] args)
		{
		File dir=new File("/Users/rodri/Desktop/goa/nuevos");
		for(File f:dir.listFiles())
			{
			String s=f.getAbsolutePath();
			if(s.contains("gene_association.goa"))
				GOparser.map(s, "/Users/rodri/desktop/goa/mini/gobiological_process-"+s.substring(s.indexOf(".goa")+5)+".map");
			}
		//GOparser.map("/Users/rodri/Documents/investigacion/distros/git/voronto/Voronto/data/calbicans/annotations/gene_association.cgd", "/Users/rodri/Documents/investigacion/distros/git/voronto/Voronto/src/es/usal/voronto/data/go/gobiological_process-calbicans_gene_ensembl.map");
		}
	
	private static List<OntologyTerm> lot;

	/**
	 * Unlike KO parser, this method only builds the hierarchy, but does not specify the gene mapping
	 * Then, method annotate should be used for gene mapping, providing the species name
	 * @param path
	 * @param ontology - either biological_process, cellular_component or molecular_function
	 * @return
	 */
	public static TreeMap<OntologyTerm, TreeMap> parse(String path, String ontology, boolean redo)
	{
	TreeMap<OntologyTerm, TreeMap> map=new TreeMap<OntologyTerm, TreeMap>();
	//Map<OntologyTerm, String> parents=Collections.synchronizedMap(new TreeMap<OntologyTerm, String>());
	Map<OntologyTerm, List<String>> parents=Collections.synchronizedMap(new TreeMap<OntologyTerm, List<String>>());//NOTE: considering more than 1 parent per node
	Map<String, List<String>> children=Collections.synchronizedMap(new TreeMap<String, List<String>>());
	
	
	String serName=null;
	if(path.contains("slim"))	serName="goSLIM";
	else						serName="go";
	
	String goTermName="terms";
	if(path.contains("slim"))	goTermName="goSLIM"+ontology+".txt";
	else						goTermName="go"+ontology+".txt";
	ArrayList<String> goTerms=new ArrayList<String>();
	
	//1) ###################################### HIERARCHY ##################################
	//A) Try to get GO hierarchy map from serialized data, if it exists
	if(!redo)
		try
		   {
			System.out.println("Reading GO serialized map");
			long time=System.currentTimeMillis();
			InputStream is=Thread.currentThread().getContextClassLoader().getResourceAsStream("es/usal/voronto/data/go/"+serName+""+ontology+".ser");
			ObjectInputStream oin = new ObjectInputStream(is);
			map = (TreeMap<OntologyTerm, TreeMap>)oin.readObject();
		    oin.close();
		    System.out.println("GO hierarchy read in "+(System.currentTimeMillis()-time)/1000.0);
		    }catch(Exception e){e.printStackTrace(); redo=true;}
	
	if(map.size()==0)
		{
		//B) Else, parse OBO file to generate hierarchy
		BufferedReader in;
		try {
			in = new BufferedReader(new InputStreamReader(Thread.currentThread().getContextClassLoader().getResourceAsStream(path)));
			
			String s=null;
			String name=null;
			String id=null;
			String namespace=null;
			String def=null;
			
			//B1) map ontology terms to their parents
			long time=System.currentTimeMillis();
			System.out.println("Building GO ontology hierarchy: "+path+"\t"+ontology);
			while ((s = in.readLine()) != null)
				{
				if(s.equals("[Term]"))	
					{
					//B11) Take first three fixed fields (id, name and namespace)
					id=in.readLine();
					id=id.substring(id.indexOf(":")+1).trim();
					name=in.readLine();
					name=name.substring(name.indexOf(":")+1).trim();
					namespace=in.readLine();
					namespace=namespace.substring(namespace.indexOf(":")+1).trim();
					OntologyTerm ot=new OntologyTerm(name, id);
					if(id.contains("0008652"))
						System.out.println("cellular aa");
					if(namespace.equals(ontology))	//only for the ontology we are populating
						{
						//B12) Then go to is_a relations
						String cad=null;
						boolean obsolete=false;
						
						//considering more than one parent
						List<String> p=new ArrayList<String>();
						while((cad=in.readLine()).length()>0)	
							{
							if(cad.startsWith("is_obsolete: true"))	//this is a non valid node
								{	obsolete=true; break;	}
							if(cad.startsWith("is_a:"))	
								p.add(cad.substring(cad.indexOf(":")+1, cad.indexOf("!")).trim());//this is the root, no is_a: relationship
							}
						
	
						//TODO: by now, we only consider one parent, as in GO slim
						if(!obsolete)
							{
							goTerms.add(id);
							parents.put(ot, p);
							
							//with children list here
							for(String pp:p)
								{
								List<String> c=children.get(pp);
								if(c!=null)	c.add(ot.id);
								else	{c=new ArrayList<String>(); c.add(ot.id); children.put(pp, c);}
								}
							}
						//considering just one parent (first one)
						/*while(!(cad=in.readLine()).startsWith("is_a:"))	
							{
							if(cad.startsWith("is_obsolete: true"))	//this is a non valid node
								{	obsolete=true; break;	}
							if(cad.length()==0)	break;//this is the root, no is_a: relationship
							}
	
						//TODO: by now, we only consider one parent, as in GO slim
						if(!obsolete)
							{
							String parentid=null;
							goTerms.add(id);
							if(cad.length()>0)
								parentid=cad.substring(cad.indexOf(":")+1, cad.indexOf("!")).trim();
							parents.put(ot, parentid);
							}*/
						}
					
	
					}//if [Term]
				}
			
			System.out.println("Ontology terms recovered in "+(System.currentTimeMillis()-time)/1000.0+"s, "+parents.size()+" elements");
			//B2) get the one without a parent (root) and build from there on
			Iterator<OntologyTerm> it=parents.keySet().iterator();
			OntologyTerm root=null;
			while(it.hasNext())
				{
				OntologyTerm ot=it.next();
				if(ot.name.equals(ontology))
				//if(parents.get(ot)==null)//da root
					{
					root=ot;
					map.put(ot, new TreeMap<OntologyTerm, Map>());
					//break;
					}
				if(parents.get(ot)==null)
					System.out.println(ot.id+" "+ot.name+" is a root");
				}
			//only one should be in the map by now
			System.out.println("Parent is "+root.id+", "+root.name);
			OntologyTerm[] ots=parents.keySet().toArray(new OntologyTerm[0]);
			ArrayList<OntologyTerm> lot=new ArrayList<OntologyTerm>(Arrays.asList(ots));
			HashMap<String, OntologyTerm> mlot=new HashMap<String, OntologyTerm>();
			for(OntologyTerm ot:ots)	mlot.put(ot.id, ot);
			
//			mlot.remove(root);
			//System.out.println("root is at "+lot.indexOf(root));
			//lot.remove(lot.indexOf(root));
			//addChildren(root, parents, map.get(root), lot,0);
			addChildren(root, children, map.get(root), mlot,0);
		//	public static synchronized void addChildren(OntologyTerm root, Map<String, List<String>> children, Map<OntologyTerm, Map> m, Map<String, OntologyTerm> lot)
				
			
				
		if(redo)
			{
			//B3) Serialize the mapping for optimization (still it takes time to read after, better divide by species
			FileOutputStream fos = new FileOutputStream(serName+ontology+".ser");
			ObjectOutputStream out = new ObjectOutputStream(fos);
			out.writeObject(map);
			out.close();
			
			//B4) Write Go terms to a file as a list:
			BufferedWriter fichero=new BufferedWriter(new FileWriter(goTermName));
			for(String go:goTerms)	fichero.write(go+"\n");
			fichero.close();
			}

				} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e){e.printStackTrace();}
	}
	
	TreeMap<OntologyTerm, TreeMap> onto=new TreeMap<OntologyTerm, TreeMap>();
	onto.put(new OntologyTerm("GO slim ontology", ""), map);
	
	return onto;
	}
	
	/**
	 * Recursively builds the hierarchy, given the parents.
	 * It cannot be done from the beginning because of the GO file structure.
	 * @param root
	 * @param parents
	 * @param m
	 * @param lot
	 */
	public static synchronized void addChildren(OntologyTerm root, Map<String, List<String>> children, Map<OntologyTerm, Map> m, Map<String, OntologyTerm> lot, int level)
		{
		//System.out.println();
		//for(int i=0;i<level;i++)	System.out.print("\t");
		//System.out.print("Adding children to "+root.name);
		//if(root.id.equals("GO:0008152"))
		//	System.out.println("cellular proc");
		List<String> c=children.get(root.id);
		if(c!=null)
			{
			for(String cc:c)
				{
				OntologyTerm ochild=lot.get(cc);
				m.put(ochild, new TreeMap<OntologyTerm, Map>());
				addChildren(ochild, children, m.get(ochild), lot, level+1);
				}
			}
		return;
		}

	/**
	 * Takes a GAF file and generatesmi a map from it (basically the same file, but smaller)
	 * 
	 * @param gafFilePath
	 */
	public static void map(String gafFilePath, String mapFilePath)
		{
		BufferedReader in;
		BufferedWriter out;
		try {
			in = new BufferedReader(new FileReader(gafFilePath));
			out = new BufferedWriter(new FileWriter(mapFilePath));
			out.write("go_id\texternal_gene_id");	out.newLine();
			String cad=null;
			while((cad=in.readLine())!=null)
				{	
				if(!cad.startsWith("!"))
					{
					String[] fields=cad.split("\t");
					String go=fields[4];
					String orf=fields[10];
					if(orf.indexOf("|")>0)	orf=orf.substring(0, orf.indexOf("|"));
					out.write(go+"\t"+orf); out.newLine();
					}
				}
			out.close();
		}catch(Exception e){e.printStackTrace();}
		}
	
	
	/**
	 * Returns the map m, but with its ontology terms mapped to gene ids on the corresponding format (entrezgene, external_gene_id, ensembl_gene_id or hgnc_symbol)
	 * This method reads annotation files produced by biomaRt via R (see script es.usal.voronto.rcode.retrieveGOannotations.R)
	 * These files are pregenerated and updated periodically (set and automate period... six months?)
	 * @param m			GO hierarchy to annotate
	 * @param species	full standard species name (e.g. Homo sapiens, but not H sapiens, hsa or human)
	 * @param gene_id	either entrezgene, external_gene_id, ensembl_gene_id or hgnc_symbol (standard id names from biomaRt)
	 * TODO: either generate a file for each GO type (Slim or not; BP/MF/CC) or big annotation files with everything and select depending if they appear in m
	 * @return
	 */
	public static TreeMap<OntologyTerm, TreeMap> annotate(TreeMap<OntologyTerm, TreeMap> m, String species, String gene_id, int type, ExpressionData ed) throws Exception 
		{
		String ontology="";
		switch(type)
			{
			case VoronoiVisualization.SLIMBP:
				ontology="goSLIMbiological_process";
				break;
			case VoronoiVisualization.SLIMCC:
				ontology="goSLIMcellular_component";
				break;
			case VoronoiVisualization.BP:
				ontology="gobiological_process";
				break;
			}
		
		BufferedReader in = null; 
		
		try{
		String path="es/usal/voronto/data/go/"+ontology+"-"+convert(species)+"_gene_ensembl.map";//by now testing with GO slim bp
		in = new BufferedReader(new InputStreamReader(Thread.currentThread().getContextClassLoader().getResourceAsStream(path)));//para applet/jws
		}catch(Exception e)
			{
			throw new Exception("No GO mapping for organism "+species);
			}
		
		String header=in.readLine();
		String []fields=header.split("\t");
		int id=-1;//column with the requested gene ids.
		for(int i=0;i<fields.length;i++)
			if(fields[i].equals(gene_id))	{id=i; break;}
		
		if(id==-1)
			{
			String s="Gene id '"+gene_id+"' is not supported for this GO map, please use one of these:\n";
			for(int i=1;i<fields.length;i++)	s+="\t"+fields[i];
			
			throw new Exception(s);
			}
		
		String cad;
		TreeMap<String, ArrayList<String>> annotations=new TreeMap<String, ArrayList<String>>();
		TreeSet<String> mappedGenes=new TreeSet<String>();
		while((cad=in.readLine())!=null)
			{
			fields=cad.split("\t");
			if(fields!=null && fields[id].length()>0 && !fields[id].equals("NA"))	//add the annotation to a data structure
				{
				if(annotations.get(fields[0])==null)
					annotations.put(fields[0], new ArrayList<String>());
				annotations.get(fields[0]).add(fields[id]);
				mappedGenes.add(fields[id].toLowerCase());
				}
			}
		
		//---
		if(ed!=null)
			{
			int nomap=0;
			
			for(String mg:ed.sortedGeneNames)
				if(!mappedGenes.contains(mg))
					{
					nomap++;
					}
			System.out.println("Genes in the expression data but not annotated: "+nomap);
			}
		//---
		
		long time=System.currentTimeMillis();
		addRecursiveAnnotations(annotations,m,null,0);
		
		System.out.println("time in adding annotations "+(System.currentTimeMillis()-time)/1000.0);
		
		return m;
		}
	
	/**
	 * Adds to each OntologyTerm in m (recursively) genes related to it, following annotations
	 * NOTE: if an OntologyTerm is annotated with some genes but its parent is not, the parent appears but without annotated genes
	 * @param annotations
	 * @param m
	 * @param level
	 */
	public static void addRecursiveAnnotations(TreeMap<String, ArrayList<String>> annotations, TreeMap<OntologyTerm, TreeMap> m, OntologyTerm parent, int level)
		{
		HashSet<String> hs=new HashSet();
		
		for(OntologyTerm ot:m.keySet())
			{
			ArrayList<String> annot=annotations.get(ot.id);
			
			if(annot!=null && ot.geneIds.size()==0)	ot.geneIds.addAll(annot);	//add annotations only once
			
			TreeMap<OntologyTerm,TreeMap> mc=m.get(ot);			//and keep recursing
			if(mc!=null)	
				addRecursiveAnnotations(annotations, mc, ot, level+1);
			
			if(parent!=null)			//and check if the parent has it (because some ontologies are inconsistent!)	NOTE: superSLOW in GO-full+mice (70s)
				parent.geneIds.addAll(ot.geneIds);
			}
		return;
		}
	
	/**
	 * Converts standard species names to OBO species names
	 * E.g. Homo sapiens --> hsapiens
	 * @param species
	 * @return
	 */
	public static String convert(String species)
		{
		return species.trim().toLowerCase().charAt(0)+species.substring(species.indexOf(" ")+1);
		}
	
	
	/**
	 * RESTful access, get annotations for a given term
	 * http://amigo.geneontology.org/cgi-bin/amigo/term-assoc.cgi?term=GO:0007200&term_assocs=direct&format=go_assoc
	 * term=GO:007200
	 * term_assocs=direct
	 * format=go_assoc
	 * 
	 * Get annotations for a given gene (GB: for entrezgene, ENSEMBL for ensembl)
	 * http://amigo.geneontology.org/cgi-bin/amigo/gp-assoc.cgi?gp=FB:FBgn0015946&format=go_assoc
	 * 
	 * So, por GO BP (full) we do have 22000 terms.
	 * We can
	 * 1) REST Search for each term its annotated genes --> 22000 queries
	 * 2) REST Search for each gene its annotated terms --> as much queries as genes in our matrix (similar to above)
	 * 3) Offline BiomaRt search for each term --> one for each organism, store large files
	 * 
	 * Pros and cons
	 * 		Memory
	 * 1	No memory increase
	 * 2	No memory increase
	 * 3	Large files storage: increase the size of the program up to 100M, even if only some organisms (human, mouse, yeast) recorded
	 * 
	 * 		Time
	 * 1	22000xt, being t the time for a REST query
	 * 2	Nxt, being N the number of genes in the expression data. Expected to be similar that in 1
	 * 3	Offline time, no charge on client
	 * 
	 * 		Update
	 * 1	Always up to date
	 * 2	Always up to date
	 * 3	Need of periodical updates
	 * 
	 * 		Memory	Time	Update
	 * 1	Good	Bad		Good
	 * 2	Good	Bad		Good
	 * 3	Bad		Good	Bad
	 * 
	 * But time is much more important in a visualizaiton, and is already burdened by tessellation.
	 * A solution like 3 is the one that has been selected for KEGG and GO slim.
	 * 
	 * Let's quantize how much it takes each thing:
	 * Time (term annotation retrieval)
	 * 1
	 * 2
	 * 3	3h15' for all the BP terms for mus musculus (4.6M) --> which means it's unafordable for running time at all
	 * 
	 *  Time (term annotation file parsing and mapping)
	 *  1
	 *  2
	 *  3
	 */
}
