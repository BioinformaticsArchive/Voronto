package es.usal.voronto.model.ontology;

import java.io.BufferedReader;
import java.io.BufferedWriter;
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
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import edu.emory.mathcs.backport.java.util.Arrays;
import es.usal.voronto.view.VoronoiVisualization;

import keggapi.Definition;
import keggapi.KEGGLocator;
import keggapi.KEGGPortType;

public class GOparser {
	
	
	
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
	Map<OntologyTerm, String> parents=Collections.synchronizedMap(new TreeMap<OntologyTerm, String>());
	
	
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
			InputStream is=Thread.currentThread().getContextClassLoader().getResourceAsStream("es/usal/voronto/data/go/"+serName+""+ontology+".ser");
			ObjectInputStream oin = new ObjectInputStream(is);
			map = (TreeMap<OntologyTerm, TreeMap>)oin.readObject();
		    oin.close();
		    System.out.println("GO hierarchy read");
		    }catch(Exception e){e.printStackTrace();}
	
	if(map.size()==0)
		{
		//B) Else, parse OBO file to generate hierarchy
		BufferedReader in;
		try {
			in = new BufferedReader(new InputStreamReader(Thread.currentThread().getContextClassLoader().getResourceAsStream(path)));//for applet/jws
			
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
					if(namespace.equals(ontology))	//only for the ontology we are populating
						{
						//B12) Then go to is_a relations
						String cad=null;
						boolean obsolete=false;
						while(!(cad=in.readLine()).startsWith("is_a:"))	
							{
							if(cad.startsWith("is_obsolete: true"))	//this is a non valid node
								{	obsolete=true; break;	}
							if(cad.length()==0)	break;//this is the root, no is_a: relationship
							}
	
						//by now, we only consider one parent, as in GO slim
						if(!obsolete)
							{
							String parentid=null;
							goTerms.add(id);
							if(cad.length()>0)
								parentid=cad.substring(cad.indexOf(":")+1, cad.indexOf("!")).trim();
							parents.put(ot, parentid);
							}
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
				if(parents.get(ot)==null)//da root
					{
					root=ot;
					map.put(ot, new TreeMap<OntologyTerm, Map>());
					}
				}
			//only one should be in the map by now
			System.out.println("Parent is "+root.id+", "+root.name);
			OntologyTerm[] ots=parents.keySet().toArray(new OntologyTerm[0]);
			ArrayList<OntologyTerm> lot=new ArrayList<OntologyTerm>(Arrays.asList(ots));
		
			System.out.println("root is at "+lot.indexOf(root));
			lot.remove(lot.indexOf(root));
			addChildren(root, parents, map.get(root), lot,0);
			
			
				
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
	public static synchronized void addChildren(OntologyTerm root, Map<OntologyTerm, String> parents, Map<OntologyTerm, Map> m, ArrayList<OntologyTerm> lot, int level)
		{
		System.out.println("Adding chidren to "+root.name);
		for(OntologyTerm ot:lot)
			{
			if(parents.get(ot).equals(root.id))
				{
				//System.out.println("Adding parent "+root.name+" to term "+ot.name);
				m.put(ot, new TreeMap<OntologyTerm, Map>());
				
				//copy list to avoid recursive concurrent modifications
				ArrayList<OntologyTerm> nlot=new ArrayList<OntologyTerm>(Arrays.asList(new OntologyTerm[lot.size()]));
				Collections.copy(nlot, lot);
				nlot.remove(nlot.indexOf(ot));
				
				addChildren(ot, parents, m.get(ot), nlot, level+1);
				}
			}
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
	public static TreeMap<OntologyTerm, TreeMap> annotate(TreeMap<OntologyTerm, TreeMap> m, String species, String gene_id, int type) throws Exception 
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
		
		String path="es/usal/voronto/data/go/"+ontology+"-"+convert(species)+"_gene_ensembl.map";//by now testing with GO slim bp
		BufferedReader in = new BufferedReader(new InputStreamReader(Thread.currentThread().getContextClassLoader().getResourceAsStream(path)));//para applet/jws
		
		String header=in.readLine();
		String []fields=header.split("\t");
		int id=-1;//column with the requested gene ids.
		for(int i=0;i<fields.length;i++)
			if(fields[i].equals(gene_id))	{id=i; break;}
		
		String cad;
		TreeMap<String, ArrayList<String>> annotations=new TreeMap<String, ArrayList<String>>();
		while((cad=in.readLine())!=null)
			{
			fields=cad.split("\t");	
			if(fields[id].length()>0 && !fields[id].equals("NA"))	//add the annotation to a data structure
				{
				if(annotations.get(fields[0])==null)
					annotations.put(fields[0], new ArrayList<String>());
				annotations.get(fields[0]).add(fields[id]);
				}
			}
		
		addRecursiveAnnotations(annotations,m,null,0);
		
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
		
		for(OntologyTerm ot:m.keySet())
			{
			if(annotations.get(ot.id)!=null)
				{
				ot.geneIds.addAll(annotations.get(ot.id));	//add annotations...
				}
			//else
			//	System.out.println(ot.name+": no annotations found");
			
			
			TreeMap<OntologyTerm,TreeMap> mc=m.get(ot);			//and keep recursing
			if(mc!=null)	addRecursiveAnnotations(annotations, mc, ot, level+1);
			
			//... and check parent's annotation to solve inconsistencies (might happen in GO!)
			//This takes a wild time, as it adds up thousands of new annotations! it also reflects in timing for voronoi final mapping
			if(parent!=null)
				{
				for(String cg:ot.geneIds)
					{
					if(!parent.geneIds.contains(cg))	
						{
						parent.geneIds.add(cg);
						}
					}
				}
			//---------------
			}
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
