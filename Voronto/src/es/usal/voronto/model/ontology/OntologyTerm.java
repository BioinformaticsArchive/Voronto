package es.usal.voronto.model.ontology;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;

public class OntologyTerm implements Serializable, Cloneable, Comparable<OntologyTerm> {
	public String name;
	public String id;
	//public ArrayList<String> geneIds=new ArrayList<String>(); //gene ids annotated with this term;
	public HashSet<String> geneIds=new HashSet<String>(); //gene ids annotated with this term;
	public Map<String, ArrayList<Float>> geneExs=new TreeMap<String, ArrayList<Float>>();//expression for each gene term, on each condition
	public Map<String, ArrayList<String>> koIds=new TreeMap<String, ArrayList<String>>();//ko ids for kegg orthology, mapped to Kegg genes
	public Map<String, ArrayList<Float>> koExs=new TreeMap<String, ArrayList<Float>>();//expression for each ko term, on each condition
	
	public OntologyTerm(String name, String id)
		{
			this.name=name;
			this.id=id;
		}

	@Override
	public int compareTo(OntologyTerm arg0) {
		// TODO Auto-generated method stub
		return 10*name.compareTo(arg0.name)+id.compareTo(arg0.id);
		}
	
}
