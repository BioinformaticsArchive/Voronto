package es.usal.voronto.model.ontology;

import java.io.BufferedReader;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import keggapi.*;
import javax.xml.rpc.ServiceException;

import es.usal.voronto.model.expression.ExpressionData;

/**
 * Parsers KEGG orthology file into a Map
 * @author rodri
 *
 */
public class KOparser {

	public static ArrayList<String> getPathways(String organism)
		{
		ArrayList<String> paths=new ArrayList<String>();
		
		try{
			KEGGLocator  locator = new KEGGLocator();
			KEGGPortType serv    = locator.getKEGGPort();
			Definition[] def=serv.list_pathways(organism);
			for(Definition d:def)
				{
				paths.add(d.getEntry_id().substring(8));
				}
			}catch(Exception e){e.printStackTrace(); return null;}
		
		return paths;
		}
	/**
	 * 
	 * 
	 * @param path KEGG hierarchy serialized file
	 * @param organism KEGG organism, so it selects only entries for that organism (e.g. if Saccharomyces cerevisiae --> sce)
	 * @param pathIDs IDs of the paths that must be included in the ontology (if null, all of them will be included)
	 * @param ed  ExpressionData which will determine selected gene annotations and gene id (it can be 'entrezgene', 'ensembl_gene_id' or 'external_gene_id')
	 * 					By default is null, and gene annotations will be the ones coming directly from KEGG (usually entrezgene, but not always)
	 * @param update If true, it uses programmatic remaps elements and updates the serialized file (takes time)
	 * 
	 * @return
	 */
	public static TreeMap<OntologyTerm, TreeMap> parse(String path, String organism, ArrayList<String> pathIDs, ExpressionData ed, boolean update)
	{
	TreeMap<OntologyTerm, TreeMap> map=new TreeMap<OntologyTerm, TreeMap>();
	TreeMap<String, ArrayList<String>> koterms=new TreeMap<String, ArrayList<String>>();//mapping of KO terms to gene ids
	
	FileInputStream fis = null;
	ObjectInputStream oin = null;
	boolean addElements=true;
	HashMap<String, String> geneMap=null;
	
	try
	   {
		System.out.println("Reading komappingTOTAL");
		InputStream is=Thread.currentThread().getContextClassLoader().getResourceAsStream("es/usal/voronto/data/kegg/komappingTOTAL.ser");
		oin = new ObjectInputStream(is);
        koterms = (TreeMap<String, ArrayList<String>>)oin.readObject();
	    oin.close();
	    System.out.println("koterms read");
	    if(ed!=null)
	    	{
	    	if(ed.organismKegg.equals("sce") || ed.organismKegg.equals("spo"))	
	    		{if(!ed.chip.equals("ensembl_gene_id") && ed.ensembl_gene_idHash!=null)	
	    			geneMap=ed.invertedHash(ed.ensembl_gene_idHash);}
	    	else								
	    		{if(!ed.chip.equals("entrezgene") && ed.entrezgeneHash!=null)	
	    			geneMap=ed.invertedHash(ed.entrezgeneHash);}
	    	}
	    }catch(Exception e){e.printStackTrace();}
	
	BufferedReader in;
	int numPaths=0;
	
	try {
		in = new BufferedReader(new InputStreamReader(Thread.currentThread().getContextClassLoader().getResourceAsStream(path)));//para applet/jws
		String s=null;
		String name=null;
		String id=null;
		OntologyTerm currentA=null, currentB=null, currentC=null, currentD=null;
		
		//NOTE: only for retrieving info (local execution,not applet)
		KEGGPortType serv=null;
		if(update)
			{
			KEGGLocator  locator = new KEGGLocator();
			serv    = locator.getKEGGPort();
			}
		
        
		while ((s = in.readLine()) != null)
			{
			switch(s.charAt(0))
				{
				case 'A':
					name=s.substring(s.indexOf("<b>")+3, s.indexOf("</b>")).trim();
					System.out.println(name);
					currentA=new OntologyTerm(name, "");
					map.put(currentA, new TreeMap<OntologyTerm, TreeMap>());
					break;
				case 'B':
					if(s.length()==1)	break;
					
					name=s.substring(s.indexOf("<b>")+3, s.indexOf("</b>")).trim();
		//			System.out.println("\t"+name);
					currentB=new OntologyTerm(name, "");
					map.get(currentA).put(currentB, new TreeMap<OntologyTerm, TreeMap>());//Cannot be done with ontologies!!!
					break;
				case 'C'://pathway
					name=s.substring(s.indexOf("0")+5).trim();
					id=s.substring(s.indexOf("0"),s.indexOf("0")+5);
					if(name.contains("["))	name=name.substring(0, name.indexOf("["));
					currentC=new OntologyTerm(name, id);
	//				System.out.println("\t\t"+name+"\t"+id);
					//TODO ignore if not the one and so on
					if(pathIDs==null || pathIDs.contains(id))	
						{
						((Map<OntologyTerm, TreeMap>)((Map<OntologyTerm, TreeMap>)map.get(currentA)).get(currentB)).put(currentC, new TreeMap<OntologyTerm, Map>());
						//System.out.println("path "+name+" available for this organism");
						numPaths++;
						addElements=true;
						}
					else
						{
						addElements=false;
						//System.err.println("path "+name+" not available for this organism");
						}
					
					break;
				case 'D'://element
					if(!addElements)	break;
					
					name=s.substring(s.indexOf("K")+6).trim();
					id=s.substring(s.indexOf("K"), s.indexOf("K")+6).trim();
					if(name.contains("["))	name=name.substring(0, name.indexOf("["));
		//			System.out.println("\t\t\t"+name);
					
					ArrayList<String> gids=new ArrayList<String>();
					if(koterms.containsKey("ko:"+id))
						gids=koterms.get("ko:"+id);	
					//NOTE: only for retrieving KO info, not for applet
					else
						{
						System.err.println("ko:"+id+" not found, should be added");
						if(update)
							{
							System.out.println("\tadding...");
							Definition[] results = serv.get_genes_by_ko("ko:"+id, "all");
							for(Definition r:results)
								gids.add(r.getEntry_id());
							System.out.println(name+": found "+results.length+" results");
							koterms.put("ko:"+id, gids);
							}
						}
					
					//TRANSLATION to gene ids on the expression data
					if(geneMap!=null)
						{
						ArrayList<String> gids2=new ArrayList<String>();
						for(String gid:gids)
							gids2.add(gid.substring(0, gid.indexOf(":"))+":"+geneMap.get(gid.substring(gid.indexOf(":")+1).toLowerCase()));
						gids=gids2;
						}
					
					
					//This should be something recursive on deeper ontologies
					Iterator<OntologyTerm> it=map.keySet().iterator();
					while(it.hasNext())
						{
						OntologyTerm t=it.next();
						if(t.name.equals(currentA.name))
							{
							if(organism==null)
								{
								t.geneIds.addAll(gids);
								t.koIds.put(id, gids);
								}
							else
								{
								ArrayList<String> l=new ArrayList<String>();
								for(String gid:gids)
									if(gid.startsWith(organism))
										{
										t.geneIds.add(gid);
										l.add(gid);
										}
								t.koIds.put(id, l);
								}
							break;
							}
						}
					
					Map<OntologyTerm, TreeMap> ma=map.get(currentA);
					it=ma.keySet().iterator();
					while(it.hasNext())
						{
						OntologyTerm t=it.next();
						if(t.name.equals(currentB.name))
							{
							if(organism==null)
								{
								t.geneIds.addAll(gids);
								t.koIds.put(id, gids);
								}
							else
								{
								ArrayList<String> l=new ArrayList<String>();
								for(String gid:gids)
									if(gid.startsWith(organism))
										{
										t.geneIds.add(gid);
										l.add(gid);
										}
								t.koIds.put(id, l);
								}
							break;
							}
						}
					Map<OntologyTerm, TreeMap> mb=ma.get(currentB);
					it=mb.keySet().iterator();
					while(it.hasNext())
						{
						OntologyTerm t=it.next();
						if(t.name.equals(currentC.name))
							{
							if(organism==null)
								{
								t.geneIds.addAll(gids);
								t.koIds.put(id, gids);
								}
							else
								{
								ArrayList<String> l=new ArrayList<String>();
								for(String gid:gids)
									if(gid.startsWith(organism))
										{
										t.geneIds.add(gid);
										l.add(gid);
										}
								t.koIds.put(id, l);
								}
							break;
							}
						}
					
					break;
				default:
					break;
				}
			}
		if(update)
			{
			System.out.println("Saving the whole thing");
			//Serialize the mapping for optimization
			FileOutputStream fos = new FileOutputStream("komappingTOTAL-new.ser");
			ObjectOutputStream out = new ObjectOutputStream(fos);
			out.writeObject(koterms);
			out.close();	
			}
		
	} catch (FileNotFoundException e) {
		e.printStackTrace();
	} catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	} catch (Exception e){e.printStackTrace();}
	TreeMap<OntologyTerm, TreeMap> onto=new TreeMap<OntologyTerm, TreeMap>();
	onto.put(new OntologyTerm("KEGG Orthology", ""), map);
	System.out.println("Number of pathways parsed: "+numPaths);
	return onto;
	}

}

