package es.usal.voronto.control;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Insets;
import java.awt.Toolkit;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Comparator;
import java.util.TreeMap;

import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.ProgressMonitor;
import javax.swing.SwingWorker;

import es.usal.voronto.model.expression.ExpressionData;
import es.usal.voronto.model.ontology.OBOparser;
import es.usal.voronto.model.ontology.ReactomeParser;
import es.usal.voronto.model.ontology.GOparser;
import es.usal.voronto.model.ontology.KOparser;
import es.usal.voronto.model.ontology.OntologyTerm;
import es.usal.voronto.model.ontology.XMLparser;
import es.usal.voronto.model.voronoi.Cell;
import es.usal.voronto.view.Introduction;
import es.usal.voronto.view.VoronoiVisualization;



/**
 * This class is the main class for Voronto, it first shows a presentation screen, then asks for input, show progress and finally
 * represent the visualization vi VoronoiVisualization.java
 * @author rodri
 *
 */


public class Voronto extends JFrame implements PropertyChangeListener{

	private static final long serialVersionUID = 5076238796400391965L;
	VoronoiVisualization gv;
	Introduction intro;
	private ProgressMonitor progressMonitor;
	private JTextArea taskOutput;
	private Task task;
    public float unit=0;
    private String message="";
    private int width=700;
    private int height=600;
    
    public String expressionPath="./.voronto";
	public String ontologyPath="./.vorontology";
	
  	
    /**
	 * Class that reads expression and ontology files, and generates tessellation
	 */
	class Task extends SwingWorker<Void,Void>
		{
		
		int ontology=0;
		//String expressionFile="";
		File expressionFile=null;
		Voronto voronto=null;
		String ontologyFile="";
		public Task(File ef, String of, int o, Voronto v)
			{
			expressionFile=ef;
			ontologyFile=of;
			ontology=o;
			voronto=v;
			}
		
		 @Override
	        public Void doInBackground() {
	            try {
			 	Long time=System.currentTimeMillis();
			 	Long time2;
			 	message=null;
			 	this.
	            setProgress(1);
	           // message="Building Voronto visualization";
	           // setProgress(5);
	            Thread.sleep(100);
	            
	           	//1) Read expression data
		        ExpressionData md=null;
		        if(expressionFile!=null)
			        {
		        	message="Loading expression data...";
			        setProgress(10);    
			        Thread.sleep(100);
			           
			        System.out.println("building expression data...");
			    	time=System.currentTimeMillis();
			    	try{
			    		md=new ExpressionData(expressionFile.getAbsolutePath(), expressionFile.getName(), false, 1,1,1);
				    	time2=System.currentTimeMillis();
				        message="\t\tDONE\tin "+(time2-time)/1000.0+"s\nParsing ontology...";
				        setProgress(20);
					    System.out.println("Time in getting the microarray: "+(time2-time)/1000.0);
					    System.out.println("Average: "+md.average+" min: "+md.min+" max: "+md.max);
		            	}
			    	catch(Exception e)
			    		{
			    		time2=System.currentTimeMillis();
			    		//message="\t\tERROR\n\t"+e.getMessage()+"\nNOTE: Error loading expression data, an unmapped ontology will be visualized\nParsing ontology...";
			    		message="\t\tERROR\n\t"+e.getMessage()+"\nNOTE: Error loading expression data";
			    		e.printStackTrace();
			    		setProgress(20);
			    		Thread.sleep(5000);
					    System.out.println("Time in getting the microarray: "+(time2-time)/1000.0);
					    return null;
				    	}
			    	}
		        else
		        	{
		        	//message="\nNOTE: No expression data selected, an unmapped ontology will be visualized\nParsing ontology...";
		        	message="\nNOTE: No expression data selected, please select one";
		        	setProgress(20);
			    	md=null;
		        	}
		        
		        
		        if(md==null)
			        {
			     	message="\n\nStop: no expression data loaded";
	                setProgress(0);
	                return null;
			        }
           
		        
		        //2) Parse ontology
	            int type=-1;
	            TreeMap<OntologyTerm, TreeMap> m=null;
	            String customName=null;
	            boolean expressionButNoAnnotation=false;
	            
	            switch(ontology)
            	    {
            	    case 0:
            	    	if(md!=null)       	    		
            	    		m=KOparser.parse("es/usal/voronto/data/kegg/ko00001.keg", md.organismKegg, KOparser.getPathways(md.organismKegg), md, false);
            	    	else            	    		
            	    		m=KOparser.parse("es/usal/voronto/data/kegg/ko00001.keg", null, null, null, false);
            	    		//m=KOparser.parse("es/usal/voronto/data/kegg/ko00001.keg", null, null, true);
            		    
            		    m=m.get(new OntologyTerm("KEGG Orthology", ""));
            		    System.out.println("Time in getting the ontology: "+(System.currentTimeMillis()-time)/1000.0);
            		    type=VoronoiVisualization.KEGG;
            		    break;
            	    case 1:
            	    	m=GOparser.parse("es/usal/voronto/data/go/gene_ontology_ext.obo", "biological_process", false);
            	    	//m=GOparser.parse("es/usal/voronto/data/go/gene_ontology_ext.obo", "biological_process", true);//for updates
            	    	m=m.get(new OntologyTerm("GO slim ontology", ""));
            	    	m=m.get(new OntologyTerm("biological_process", "GO:0008150"));
            	    	if(md!=null)	GOparser.annotate(m, md.organism, md.chip, VoronoiVisualization.BP, md);
            			
            	    	type=VoronoiVisualization.BP;
            			break;
            	    case 2:
            	    	m=GOparser.parse("es/usal/voronto/data/go/gene_ontology_ext.obo", "cellular_component", false);
            	    	m=m.get(new OntologyTerm("GO slim ontology", ""));
            	    	m=m.get(new OntologyTerm("cellular_component", "GO:0005575"));
            	    	if(md!=null)	GOparser.annotate(m, md.organism, md.chip, VoronoiVisualization.CC, md);
            			
            	    	type=VoronoiVisualization.CC;
            			break;	
            	    case 3:
            	    	m=GOparser.parse("es/usal/voronto/data/go/gene_ontology_ext.obo", "molecular_function", false);
            	    	m=m.get(new OntologyTerm("GO slim ontology", ""));
            	    	m=m.get(new OntologyTerm("molecular_function", "GO:0003674"));
            	    	if(md!=null)	GOparser.annotate(m, md.organism, md.chip, VoronoiVisualization.MF, md);
            			
            	    	type=VoronoiVisualization.MF;
            			break;	
            		case 4:
            			if(md!=null)  			m=ReactomeParser.readSer("es/usal/voronto/data/reactome/"+md.organism+".ser", md);
            			else       				m=ReactomeParser.readSer("es/usal/voronto/data/reactome/Homo sapiens.ser", null);
            		    type=VoronoiVisualization.REACTOME;
            		    break;
            		case 5:
            			if(intro.ontologyType.contains("xml"))	//XML custom for ontology and mapping
            				{
            				XMLparser xp=new XMLparser();
            				m=xp.parse(ontologyFile);
                			customName=xp.root.name;
                		    }
            			else		//OBO for ontology only (add GAF?)
            				{
            				//1) Parse OBO ontology + report errors
            				m=OBOparser.parse(ontologyFile, null, null);
            				//errors
            				//2) If no errors, ask for a GAF file (link to format, tell about mandatory fields)
            				if(m!=null && m.size()>0 && md!=null)
            					{
            					if(intro.annotationFile!=null && intro.annotationFile.length()>0)
                					{
            						//parse GAF file and annotate m
                    				m=GOparser.annotateFromGAF(m, intro.annotationFile.getAbsolutePath(), md);	
                					}
                				else if(md!=null)	expressionButNoAnnotation=true;
            					}
            				
            				}
            			
            			type=VoronoiVisualization.CUSTOM;
            			if(m.size()==1)	
        					m=m.get(m.keySet().iterator().next());
        				
            			break;
            	    default:
                		System.out.println("Ontology not supported");
                		break;
            	    }
                time2=System.currentTimeMillis();
                message="\t\tDONE\tin "+(time2-time)/1000.0+"s";
                if(ontology>0 && ontology<4)
                	message+="\nNOTE: Full GO is a very large ontology, only the two first levels will be visualized.\n\tIf you want to explore deeper, hover over a term and press 'enter'.";
                if(expressionButNoAnnotation)
                	message+="\nNOTE: No annotations provided for custom ontology, expression data will be ignored\n";
                	
                message+="\nBuilding tessellation...";
				setProgress(65);
        		
	            setProgress(60);
			        	
                //3) Build voronoi tessellation
                if(m!=null)
                	{
                	int maxDepth=3;			//KEGG
                	if(ontology==1 || ontology==2 || ontology==3)
                        		maxDepth=2;	//GO ontologies now only curated, no IAE (much smaller)
                	else if(ontology==4 || ontology==5)	maxDepth=20;
                	else		maxDepth=20;
                		
                	time=System.currentTimeMillis();
                	gv=null;
                	gv=new VoronoiVisualization(m, customName, md, voronto.getWidth(), voronto.getHeight(), type, maxDepth, voronto);
                	if(gv!=null)
	                	{
	            	    time2=System.currentTimeMillis();
		                message="\t\tDONE\tin "+(time2-time)/1000.0+"s\nPreparing visualization...";
		                setProgress(100);
	                	}
                	else
                		{
                		message="";
		                setProgress(0);
	                	}
	                }
                else
                	{
                	System.err.println("No hierarchy!");
                	}
            } catch (Exception e) 
            	{
            	message="\n\n ERROR: "+e.getMessage();
            	//message="\n\n ERROR: "+getError(e);
            	e.printStackTrace();
            	setProgress(0);
            	}
            return null;
	        }
	 
	        @Override
	        public void done() {
	        	progressMonitor.setProgress(0);
	        }
		}
	
	public String getError(Exception e)
		{
		String msg=e.getMessage();
		for(StackTraceElement ste:e.getStackTrace())
			msg+="\n\tat "+ste.getClassName()+"."+ste.getMethodName()+":"+ste.getLineNumber();
		return msg;
		}
	
	public static void main (String[] args) {
        
		 try{
		Voronto t = new Voronto();    // Create applet/application
		t.setSize(t.width, t.height);    // Set window size
        t.setPreferredSize(new Dimension(t.width, t.height));    // Set window size
	    t.setLayout(new BorderLayout());
	    t.setTitle("Voronto");           // Set window title
	    t.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        
	    //0) Present an empty visualization
	    t.intro=new Introduction(t);
	    t.taskOutput = new JTextArea(8, 20);
        t.taskOutput.setMargin(new Insets(5,5,5,5));
        t.taskOutput.setEditable(false);
 
        t.add(t.intro);
	    t.add(new JScrollPane(t.taskOutput), BorderLayout.SOUTH);
        t.setVisible(true);
	    t.repaint();
	   	}catch(Exception e){e.printStackTrace(); return;}
	    
    }
	
//public void launch(String expressionFile, String ontologyFile, int ontology)
public void launch(File expressionFile, String ontologyFile, int ontology)
	{
	try{
	progressMonitor = new ProgressMonitor(Voronto.this,  "Loading Voronto", "", 0, 100);
	message="";
	progressMonitor.setProgress(0);
	progressMonitor.setMillisToDecideToPopup(1000000000);
	progressMonitor.setMillisToPopup(1000000000);
	
	if(expressionPath!=null && expressionFile!=null)
		{
		BufferedWriter bw=new BufferedWriter(new FileWriter(expressionPath));
		bw.write(expressionFile.getAbsolutePath());
		bw.close();
		}
	if(ontologyPath!=null && ontologyFile!=null)
		{
		BufferedWriter bw=new BufferedWriter(new FileWriter(ontologyPath));
		bw.write(ontologyFile);
		bw.close();
		}
	
	gv=null;
	task = new Task(expressionFile, ontologyFile, ontology, this);
	task.addPropertyChangeListener(this);
	task.execute();
	
	setResizable(false);
	String path=expressionFile.getAbsolutePath();
	if(expressionFile!=null)	setTitle(path.substring(path.lastIndexOf("/")+1));
	}catch(Exception e){e.printStackTrace();}
	}
  
/**
 * Resets current voronoi visualization and switches back to introduction
 */
public void goBack()
	{
	remove(gv);
	gv=null;
	System.gc();
	
	intro.reset();
	
	add(intro);
	this.setResizable(true);
    this.setTitle("Voronto");
	
    intro.requestFocusInWindow();
    repaint();
    }

/**
 * Invoked when task's progress property changes.
 */
public void propertyChange(PropertyChangeEvent evt) {
	if ("progress" == evt.getPropertyName() ) {
        int progress = (Integer) evt.getNewValue();
        
        if(message==null)		        		taskOutput.setText("");
        else if(message.startsWith("ERROR"))	taskOutput.setText(message);
        else						        	taskOutput.append(message);
    	progressMonitor.setProgress(progress);
        if (progressMonitor.isCanceled() || task.isDone()) 
        	{
            Toolkit.getDefaultToolkit().beep();
            if (progressMonitor.isCanceled()) 
            	{
                task.cancel(true);
                taskOutput.append("Task cancelled.\n");
            	}
            else if(gv!=null)
            	{
                taskOutput.setText("");
                this.remove(intro);
                this.remove(taskOutput);
                this.add(gv);
                gv.requestFocusInWindow();
                gv.init();
            	}	
            return;
            }
    	}
	}


public Voronto(VoronoiVisualization gv)
	{
	super();
	this.gv=gv;
	}

public Voronto()
	{
	super();
	
	}


public static final class WeightComparator implements Comparator<Cell>
	{
	public int compare(Cell a, Cell b)
		{
		if(a.numLeaves>b.numLeaves)	return 1;
		else						return -1;
		}
	}
}