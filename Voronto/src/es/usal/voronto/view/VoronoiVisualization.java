package es.usal.voronto.view;

import java.awt.Color;
import java.awt.Cursor;
import java.awt.Rectangle;
import java.awt.event.KeyEvent;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeMap;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.filechooser.FileFilter;

import keggapi.KEGGLocator;
import keggapi.KEGGPortType;
import keggapi.PathwayElement;

import edu.emory.mathcs.backport.java.util.Arrays;
import es.usal.voronto.control.Voronto;
import es.usal.voronto.model.expression.ExpressionData;
import es.usal.voronto.model.ontology.OntologyTerm;
import es.usal.voronto.model.voronoi.Cell;
import es.usal.voronto.model.voronoi.BernhardtTessellation;

import processing.core.*;

public class VoronoiVisualization extends PApplet{

	/**
	 * 
	 */
	private static final long serialVersionUID = 2515531690539870969L;
	BernhardtTessellation v;
	int hoveredRegion=-1;
	int selected=-1;
	int levelThreshold=3;
	int selectedCol=0;
	private List<Cell> hoveredCells=Collections.synchronizedList(new ArrayList<Cell>());
	private List<Cell> searchedCells=Collections.synchronizedList(new ArrayList<Cell>());
	private int minHoveredLevel;
	public String message=null;
	
	public double minExp=0;
	public double maxExp=0;
	public double avgExp=0;
	public ExpressionData expData;//TODO: key controls at Voronto.java
	//scale constants to get minExp, maxExp and avgExp
	public static final int SCALE_MATRIX=0;//from the whole expression matrix
	public static final int SCALE_CONDITION=1;//from the current column
	public static final int SCALE_ONTOLOGY=2;//from the elements of the matrix in the ontology 
	public int SCALE_MODE=SCALE_CONDITION;
	
	//coloring constants for
	public static final int COLOR_EXPRESSION=0;//max or min
	public static final int COLOR_DEVIATION=1;//overall deviation respect to average
	public static final int COLOR_INTERNAL_DEVIATION=2;//overall deviation respect to average (of each term)
	public int COLOR_MODE=COLOR_EXPRESSION;
	private boolean SHOW_LABELS=true;
	
	private int contAvgExp=0;
	
	
	private int START_Y=20;//Beginning and end of Voronoi tessellation
	private int END_Y=500+START_Y;
	private Cell[] cells;
	//Boxes for each element
	Rectangle2D.Float conditionBox;
	Rectangle2D.Float organismBox;
	Rectangle2D.Float ontologyBox;
	Rectangle2D.Float scaleBox;
	Rectangle2D.Float cellBox;
	Rectangle2D.Float noteBox;
	
	private int hoveredBox=-1;
	private final int ORGANISM=0;
	private final int ONTOLOGY=1;
	private final int CONDITION=2;
	private final int CELL=3;
	private final int SCALE=4;
	public static final int KEGG=5;
	public static final int GO=6;
	public static final int SLIMBP=7;
	public static final int SLIMCC=8;
	public static final int BP=9;
	public static final int CC=10;
	public static final int REACTOME=11;
	public static final int CUSTOM=12;
	public static final int MEDIAN=0;
	public static final int MEAN=1;
	//public int whiteValue=MEDIAN;
	public int whiteValue=MEAN;
	
	private int type;
	private boolean computedLabels=false;
	
	private TreeMap<OntologyTerm,BernhardtTessellation> tessellations;//relation of tessellations we do have, in the case of redefining the tessellation by getting into it
	private OntologyTerm root; //terms that makes the root of the current tessellation
	private TreeMap<OntologyTerm, TreeMap> map;
	private int maxDepth;
	private boolean drawSearch=false;
//	private String searchString="";
	public float[] maxSd=null;
	private String helpMessage="press 'h' for help";
	private String noteLabel=helpMessage;
	private Voronto voronto=null;
	private boolean saving;
	private double medianExp;

	private String ontologyName;
	public String customOntologyName;
	
	public PFont font;
	private boolean ctrlDown=false;
	private boolean altDown=false;
	private JFrame searchFrame=null;
	private SearchFrame searchPApplet;
	private Cell minHoveredCell;
	private JFrame heatframe;
	private CellHeatmap gh;
	private int textSize=20;
	private DifExpSearch deSearchFrame;
//	private final int MAX_TEXT_SIZE=8;
//	private final int MIN_TEXT_SIZE=20;
	
	public VoronoiVisualization()
	 	{
		 super();
	 	}
	 
	public void setup()
		{
		smooth();
		noLoop();
	 	}
	
	public VoronoiVisualization(Cell[] c, int width, int height)
		{
		super();
		v=new BernhardtTessellation(c, width, height);//, this);
		}
	
	public VoronoiVisualization(TreeMap<OntologyTerm,TreeMap> m) throws Exception
		{
		super();
	 	size(width, height);
		v=new BernhardtTessellation(m, null, 700, 500,/* this,*/ VoronoiVisualization.KEGG, 3);
		cells=v.getCells();
		for(Cell c:cells)	recursiveRegionTranslation(c, 0, START_Y);
		}
	
	public VoronoiVisualization(TreeMap<OntologyTerm,TreeMap> m, String customName, ExpressionData md, int width, int height, int type, int maxDepth, Voronto parent) throws Exception
	{
	super();
	this.type=type;
	this.width=width;
	this.height=height;
	this.setSize(width, height);
	this.END_Y=height-100+START_Y;
	this.maxDepth=maxDepth;
	this.voronto=parent;
	if(customName!=null)	customOntologyName=customName;
	
 	v=new BernhardtTessellation(m, md, width, height-100, type, maxDepth);
 	this.expData=md;
 	if(type==VoronoiVisualization.KEGG)	this.levelThreshold=v.maxLevel+1;
 	if(type==VoronoiVisualization.SLIMBP || type==VoronoiVisualization.SLIMCC)	this.levelThreshold=v.maxLevel;
	cells=v.getCells();
	for(Cell c:cells)	recursiveRegionTranslation(c, 0, START_Y);
	
	tessellations=new TreeMap<OntologyTerm, BernhardtTessellation>();
	root=new OntologyTerm("root", "root");
	tessellations.put(root, v);
	map=m;
	
	setOntologyName();
	font=loadFont("es/usal/voronto/font/AppleSymbols-14.vlw");
	
	if(md!=null)	
		{
		mapExpression(md);
		expression2color(md);
		System.out.println("Finished mapping colors");
		}
	else
		expression2color(null);
	}
	
	public void setOntologyName()
		{
		if(root.name.equals("root") && root.id.equals("root"))
			{
			if(type==KEGG)	ontologyName="KEGG Orthology ("+this.levelThreshold+")";
			else if(type==REACTOME)	ontologyName="REACTOME ("+this.levelThreshold+")";
			else if(type==SLIMBP)	ontologyName="GO Slim - BP ("+this.levelThreshold+")";
			else if(type==SLIMCC)	ontologyName="GO Slim - CC ("+this.levelThreshold+")";
			else if(type==GO)	ontologyName="GO - BP ("+this.levelThreshold+")";
			else			
				{
				if(customOntologyName!=null && customOntologyName.length()>0)
					ontologyName=customOntologyName+" ("+this.levelThreshold+")";
				else
					ontologyName="Ontology ("+this.levelThreshold+")";
				}
			}
		else
			ontologyName=root.name+" ("+this.levelThreshold+")";
		}
	
	public VoronoiVisualization(int width, int height)
		{
		super();
		}
	
	public void setOntology(TreeMap<OntologyTerm, TreeMap> m) throws Exception
		{
		v=new BernhardtTessellation(m, null, width, height, /*this, */VoronoiVisualization.KEGG, 3);
		}
	
	
	/**
	 * Returns the width in pixels of a text string (cad), drawn with the given font (f)
	 * @param cad
	 * @param f
	 * @return
	 */
	public float getWidth(String cad, PFont f)
		{
		float width=0;
		for(int i=0;i<cad.length();i++)
			{
			width+=f.width(cad.charAt(i));
			}
		return width*f.getSize();
		}

	public void draw() 
		{
		fill(255,255,255);
		rect(0,0,width,height);
		noFill();
		
		if(v==null)
			{
			fill(0);
			textAlign(CENTER, CENTER);
			text(message, (int)(width*0.5), (int)(height*0.5));
			
			ellipse(50,50,100,100);
			return;
			}
		
		//drawing labels
		textFont(font, 14);
		
		fill(0,0,0);
		if(expData!=null)
			{
			textFont(font);
			textAlign(LEFT, CENTER);
			text(expData.organism, 10, 8);
			
			organismBox=new Rectangle2D.Float(10,8,getWidth(expData.organism, font), font.getSize());
			textAlign(RIGHT, CENTER);
			text(expData.getColumnLabel(selectedCol), width-10, 8);
			conditionBox=new Rectangle2D.Float(width-10-getWidth(expData.getColumnLabel(selectedCol), font),8, getWidth(expData.getColumnLabel(selectedCol), font), font.getSize());
			
			drawScale();
			//drawSkewedScale();
			}
		
		//Ontology label
		fill(0,0,0);
		textAlign(CENTER, CENTER);
		text(ontologyName, (int)(width*0.5), 8);	//TODO: by now, it's fixed!
		ontologyBox=new Rectangle2D.Float((int)(width*0.5-getWidth(ontologyName, font)*0.5),8, getWidth(ontologyName, font), font.getSize());
		
		//Note label
		if(!saving)
			{
			textAlign(CENTER, CENTER);
			fill(200,200,200);
			text(noteLabel, (int)(width*0.5-noteLabel.length()*0.5), this.END_Y+45);
			noFill();
			}
		
		textFont(font,14);
		textAlign(LEFT, CENTER);
	
		//drawing regions
		if(!computedLabels)	
			for(Cell c:cells)	
				{	
				recursiveLabelComputing(c,0);	
				computedLabels=true;
				}
		for(Cell c:cells)	recursiveRegionDrawing(c, 0);
		
		//drawing hovered region
		synchronized(hoveredCells)
			{
			for(Cell h:hoveredCells)
				{
				if(h.level<=levelThreshold)
					{
					if(h.level!=min(levelThreshold, minHoveredLevel))
						{
						stroke(0,0,0);
						strokeWeight(getWidth(h)+1);
						
						noFill();
						h.region.draw(this);
						strokeWeight(0);
						}
					else
						{
						noFill();
						strokeWeight(getWidth(h)+1);
						this.stroke(0,0,0);
						h.region.draw(this);
						fill(0,0,0);
						strokeWeight(0);
						
						textFont(font,14);
						//NOTE: numLeaves tells the number of genes in the whole of the subterms (with repetitions)
						//		h.term.geneExs.size() will give the number without repetitions.
						text(h.term.name+" ("+(int)(h.numLeaves)+")", 10, (int)(END_Y+font.getSize()*0.5));
						//text(h.term.name+" ("+h.term.geneExs.size()+")", 10, (int)(END_Y+font.getSize()*0.5));
						//text(h.term.name+" ("+h.term.geneExs.size()+" "+(int)h.numLeaves+")", 10, (int)(END_Y+font.getSize()*0.5));
						fill(255,255,255);
						
						drawProfile(h);
						}
					}
				//System.out.println("ids on term "+h.term.name+": "+h.term.geneIds.size());
				}
			if(hoveredCells.size()==0)
				{
				textFont(font,14);
				fill(0,0,0);
				text("No term selected", 10, (int)(END_Y+font.getSize()*0.5));
				noFill();
				cellBox=new Rectangle2D.Float(10,this.END_Y+8, getWidth("No term selected", font), font.getSize());
				}
			}
		
		//drawing searched regions
		for(Cell s:searchedCells)
			{
			if(COLOR_MODE==COLOR_EXPRESSION)	stroke(0,200,0);
			else								stroke(200,0,0);
			strokeWeight(getWidth(s)+3);
			
			noFill();
			if(s.region!=null)	
				{
				if(s.level<=this.levelThreshold)
					{
					boolean draw=true;
					if(s.subcells!=null)
						for(Cell s2:s.subcells)
							{
							if(searchedCells.contains(s2) && levelThreshold>=s.level+1)
								{
								draw=false;
								break;
								}
							}
						
					if(draw)	s.region.draw(this);
					}
				}
			strokeWeight(0);
			}
		
		//drawing help messages
		switch(hoveredBox)
			{
			case ORGANISM:
				drawNote("Organism of the expression data", organismBox, font);
				break;
			case ONTOLOGY:
				drawNote("Ontology in the tessellation.\nEach cell represents an ontology term.\nCell size approximates to its number of genes\nHigh level terms have wider edges.\nChange depth level with up/down arrow keys.\nThe (number) is the max depth visualized", ontologyBox, font);
				break;
			case CONDITION:
				drawNote("Expression for this condition is color-mapped.\nBrowse conditions with left/right arrow keys", conditionBox, font);
				break;
			case CELL:
				if(type==KEGG)	drawNote("Term names appear here when hovered.\nIn parenthesis, the #Êof genes in the term.\nCells color is the mean gene expression level.\nA cell profile bar is shown close to this name.\nPress 'l' to show/hide fitting labels\nDouble-click on a cell to open the browser\nwith its colored KEGG pathway.\nAlt-click on a cell to show its heatmap.", cellBox, font);
				else			drawNote("Term names appear here when hovered.\nIn parenthesis, the #Êof genes in the term.\nCells color is the mean gene expression level.\nA cell profile bar is shown close to this name.\nPress 'l' to show/hide fitting labels\nDouble-click on a cell to open the browser\nwith its ontology web entry.\nAlt-click on a cell to show its heatmap.", cellBox, font);
				break;
			case SCALE:
				if(COLOR_MODE==COLOR_EXPRESSION)	
					{
					if(whiteValue==MEDIAN)	drawNote("Quantile color scale\nPress 'm' to switch to numerical scale\n", scaleBox, font);
					if(whiteValue==MEAN)	drawNote("Numerical color scale\nPress 'm' to switch to quantile scale\n", scaleBox, font);
					
					}
				else						drawNote("Color scale for the mapping:\n  -Bright green/white for max to no deviation.\nUse 'c' key to change the scale\n", scaleBox, font);
				break;
			}
		
		//TODO: search dialog, etc. (better on a JPanel)
		//if(drawSearch)
			//text("search for: "+searchString);
		
		/*if(saving)	//Frame for pictures
			{
			strokeWeight(1);
			stroke(0,0,0);
			noFill();
			rect(0,0,width,height);
			}*/
		}

public void drawNote(String text, Rectangle2D.Float bounds, PFont font)
	{
	int margin=5;
	int noteWidth=250;
	int noteHeight=font.getSize();
	String temp=text;
	while(temp.indexOf("\n")!=-1)	
		{
		noteHeight+=font.getSize()*1.2+margin; 
		temp=temp.substring(temp.indexOf("\n")+1);
		}
	int noteX=(int)bounds.getX();
	if(bounds.x+noteWidth+margin>this.width) noteX=this.width-noteWidth-margin;	
	int noteY=(int)bounds.getY();
	if(bounds.y+noteHeight+margin+font.getSize()>this.height) noteY=this.height-noteHeight-margin-font.getSize();	
	
	fill(250,250,200);
	stroke(100,100,100);
	rect(noteX, noteY+font.getSize(), noteWidth, noteHeight);
	fill(100,100,100);
	textFont(font);
	textAlign(LEFT);
	text(text, noteX+margin, (int)(noteY+font.getSize()*1.5+margin));
	noFill();
	noStroke();
	return;
	}

/**
 * Draws profile for the given cell
 */
public void drawProfile(Cell c)
	{
	if(expData==null)	return;
	
	int squareSize=12;
	
	strokeWeight(1);
	stroke(150);
	for(int i=0;i<expData.getNumConditions();i++)
		{
		if(selectedCol!=i)
			{
			Color co=c.color.get(i);
			fill(co.getRed(), co.getGreen(), co.getBlue());
			//System.out.println(co.getRed()+" "+co.getGreen()+" "+co.getBlue());
			rect((int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+i*squareSize, this.END_Y+20,squareSize,squareSize);
			}
		}
	
	//strokeWeight(2);
	//stroke(0);
	Color co=c.color.get(selectedCol);
	fill(co.getRed(), co.getGreen(), co.getBlue());
	rect((int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+selectedCol*squareSize, this.END_Y+20,squareSize,squareSize);
	strokeWeight(2);
	stroke(0);
	line((int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+selectedCol*squareSize, this.END_Y+20+squareSize,(int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+selectedCol*squareSize+squareSize,this.END_Y+20+squareSize);
	}

/**
 * Draws a scale bar that is skewed as data, to avoid confusions
 */
public void drawSkewedScale()
	{
	int sw=150;
	int sh=10;
	int x=width-sw-15;
	int y=(int)(END_Y+sh*0.5);
	scaleBox=new Rectangle2D.Float(x, y, sw, sh);
	strokeWeight(1);
	
	Cell hc=null;
	Cell finalhc=null;
	for(int i=0;i<hoveredCells.size();i++)
		{
		hc=hoveredCells.get(i);
		if(hc.level==min(levelThreshold, minHoveredLevel))
			finalhc=hc;
		}
	hc=finalhc;
	
	boolean refDrawn=false;
	double jump0=(255-(0/(float)sw)*255);
	double jump1=(255-(0/(float)sw)*255);
	double jump=Math.ceil(jump1-jump0)+2;
	
	
	double mid=sw*0.5;//for median
	if(whiteValue==MEAN)	
		mid=sw*(expData.average-expData.min)/(expData.max-expData.min);
	
	double jumpBelow=Math.ceil(255/(255*mid/sw))+0;
	double jumpAbove=Math.ceil(255/Math.abs(255*mid/sw-255))+0;
			
	for(int i=0;i<sw;i++)
		{
		int h=0;
		boolean drawReference=false;
		switch(COLOR_MODE)
			{
			case COLOR_EXPRESSION:
				if(i>mid)
					{
					h=(int)Math.round(255-((i-mid)/mid*2*255));
					
					stroke(255,h,h);
					if(hc!=null && hc.color.get(selectedCol).getRed()==255 && Math.abs(hc.color.get(selectedCol).getBlue()-h)<=jumpAbove)
						drawReference=true;
					}
				else	
					{
					h=(int)Math.round(255-(Math.abs(i-mid)/Math.abs(mid))*255);
					stroke(h,h,255);
					if(hc!=null && hc.color.get(selectedCol).getBlue()==255 && Math.abs(hc.color.get(selectedCol).getRed()-h)<=jumpBelow)
						drawReference=true;
					}
				break;
			case COLOR_DEVIATION:
				h=(int)Math.round(255-(i/(float)sw)*255);
				stroke(h,255,h);
				if(hc!=null && Math.abs(hc.color.get(selectedCol).getBlue()-h)<=jump)
					drawReference=true;
				break;
			}
		line(x+i,y,x+i,y+sh);
		if(drawReference && !refDrawn)
			{
			stroke(0);
			line(x+i, y-2, x+i, y+sh+2);
			refDrawn=true;
			}
		}
	stroke(200,200,200);
	noFill();
	rect(x-1,y-1,sw+1,sh+1);
	
	//Draw legend
	textAlign(CENTER, TOP);
	fill(200);
	line(x,y+sh,x,y+sh+2);
	if(whiteValue==MEAN)	text("min", x, y+sh+3);
	else if(whiteValue==MEDIAN)	text("0", x, y+sh+3);
	line((int)(x+mid),y+sh,(int)(x+mid),y+sh+2);
	if(whiteValue==MEAN)		text("avg", (int)(x+mid), y+sh+3);
	else if(whiteValue==MEDIAN)	text("50", (int)(x+mid), y+sh+3);
	
	line(x+sw-1,y+sh,x+sw-1,y+sh+2);
	if(whiteValue==MEAN)		text("max", x+sw-1, y+sh+3);
	else if(whiteValue==MEDIAN)	text("100", x+sw-1, y+sh+3);
	
	line((int)(x+mid*0.5),y+sh,(int)(x+mid*0.5),y+sh+1);
	line((int)(x+mid*2),y+sh,(int)(x+mid*2),y+sh+1);
	/*line((int)(x+sw*0.1),y+sh,(int)(x+sw*0.1),y+sh+1);
	line((int)(x+sw*0.2),y+sh,(int)(x+sw*0.2),y+sh+1);
	line((int)(x+sw*0.3),y+sh,(int)(x+sw*0.3),y+sh+1);
	line((int)(x+sw*0.4),y+sh,(int)(x+sw*0.4),y+sh+1);
	line((int)(x+sw*0.6),y+sh,(int)(x+sw*0.6),y+sh+1);
	line((int)(x+sw*0.7),y+sh,(int)(x+sw*0.7),y+sh+1);
	line((int)(x+sw*0.8),y+sh,(int)(x+sw*0.8),y+sh+1);
	line((int)(x+sw*0.9),y+sh,(int)(x+sw*0.9),y+sh+1);*/
	noStroke();
	}
/**
 * Draws a scale bar
 */
public void drawScale()
	{
	int sw=150;
	int sh=10;
	int x=width-sw-15;
	int y=(int)(END_Y+sh*0.5);
	scaleBox=new Rectangle2D.Float(x, y, sw, sh);
	strokeWeight(1);
	
	Cell hc=null;
	Cell finalhc=null;
	for(int i=0;i<hoveredCells.size();i++)
		{
		hc=hoveredCells.get(i);
		if(hc.level==min(levelThreshold, minHoveredLevel))
			finalhc=hc;
		}
	hc=finalhc;
	
	boolean refDrawn=false;
	double jump0=(255-(0/(float)sw)*255);
	double jump1=(255-(0/(float)sw)*255);
	double jump=Math.ceil(jump1-jump0)+2;
	
	for(int i=0;i<sw;i++)
		{
		int h=0;
		boolean drawReference=false;
		switch(COLOR_MODE)
			{
			case COLOR_EXPRESSION:
				if(i>sw*0.5)
					{
					h=(int)Math.round(255-((i-sw*0.5)/(sw*0.5)*255));
					stroke(255,h,h);
					if(hc!=null && hc.color.get(selectedCol).getRed()==255 && Math.abs(hc.color.get(selectedCol).getBlue()-h)<=jump)
						drawReference=true;
					}
				else	
					{
					h=(int)Math.round(255-(Math.abs(i-sw*0.5)/Math.abs(sw*0.5))*255);
					stroke(h,h,255);
					if(hc!=null && hc.color.get(selectedCol).getBlue()==255 && Math.abs(hc.color.get(selectedCol).getRed()-h)<=jump)
						drawReference=true;
					}
				break;
			case COLOR_DEVIATION:
				h=(int)Math.round(255-(i/(float)sw)*255);
				stroke(h,255,h);
				if(hc!=null && Math.abs(hc.color.get(selectedCol).getBlue()-h)<=jump)
					drawReference=true;
				break;
			}
		line(x+i,y,x+i,y+sh);
		if(drawReference && !refDrawn)
			{
			stroke(0);
			line(x+i, y-2, x+i, y+sh+2);
			refDrawn=true;
			}
		}
	stroke(200,200,200);
	noFill();
	rect(x-1,y-1,sw+1,sh+1);
	
	//Draw legend
	textAlign(CENTER, TOP);
	fill(200);
	line(x,y+sh,x,y+sh+2);
	if(whiteValue==MEAN)	text("min", x, y+sh+3);
	else if(whiteValue==MEDIAN)	text("0", x, y+sh+3);
	line((int)(x+sw*0.5),y+sh,(int)(x+sw*0.5),y+sh+2);
	if(whiteValue==MEAN)		text("avg", (int)(x+sw*0.5), y+sh+3);
	else if(whiteValue==MEDIAN)	text("50", (int)(x+sw*0.5), y+sh+3);
	
	line(x+sw-1,y+sh,x+sw-1,y+sh+2);
	if(whiteValue==MEAN)		text("max", x+sw-1, y+sh+3);
	else if(whiteValue==MEDIAN)	text("100", x+sw-1, y+sh+3);
	
	//line((int)(x+sw*0.25),y+sh,(int)(x+sw*0.25),y+sh+1);
	//line((int)(x+sw*0.75),y+sh,(int)(x+sw*0.75),y+sh+1);
	line((int)(x+sw*0.1),y+sh,(int)(x+sw*0.1),y+sh+1);
	line((int)(x+sw*0.2),y+sh,(int)(x+sw*0.2),y+sh+1);
	line((int)(x+sw*0.3),y+sh,(int)(x+sw*0.3),y+sh+1);
	line((int)(x+sw*0.4),y+sh,(int)(x+sw*0.4),y+sh+1);
	line((int)(x+sw*0.6),y+sh,(int)(x+sw*0.6),y+sh+1);
	line((int)(x+sw*0.7),y+sh,(int)(x+sw*0.7),y+sh+1);
	line((int)(x+sw*0.8),y+sh,(int)(x+sw*0.8),y+sh+1);
	line((int)(x+sw*0.9),y+sh,(int)(x+sw*0.9),y+sh+1);
	noStroke();
	}


/**
 * Recursive browsing of the ontology hierarchy, drawing each cell following the tessellation placement and expression level
 * @param cell
 * @param level
 */
public void recursiveRegionDrawing(Cell cell, int level)
	{
	if(cell.region!=null && level<levelThreshold)
		{
		//1) Draw region
		this.stroke(120, 120, 120);
		strokeWeight(getWidth(cell));
		
		fill(cell.color.get(selectedCol).getRed(),cell.color.get(selectedCol).getGreen(),cell.color.get(selectedCol).getBlue());
		cell.region.draw(this); // draw this shape
		noFill();
		strokeWeight(0);
		
		if(SHOW_LABELS && cell.labelSize>0)
			{
			fill(0,0,0);
			PFont font=loadFont("es/usal/voronto/font/AppleSymbols-"+cell.labelSize+".vlw");
			
			textFont(font, cell.labelSize);
			
			text(cell.term.name, cell.labelX, cell.labelY);
			
			noFill();
			textFont(font, 14);
			}
		
		//2) Continue with recursion
		if(cell.subcells!=null && cell.subcells.length>0)	
			for(Cell cc:cell.subcells)
				recursiveRegionDrawing(cc, level+1);
		}
	}

public int getWidth(Cell cell)
	{
	//return Math.min(0, this.maxDepth*2-cell.level*2);
	return Math.max(0, 6-cell.level*2);
	}

/**
* Recursive browsing of the ontology hierarchy, computing label positions for each cell
* @param cell
* @param level
*/
public void recursiveLabelComputing(Cell cell, int level)
	{
	if(cell.region!=null)
		{
		computeLabelPosition(cell);
			
		//2) Continue with recursion
		if(cell.subcells!=null && cell.subcells.length>0)	
			for(Cell cc:cell.subcells)
				recursiveLabelComputing(cc, level+1);
			}
	}

public void computeLabelPosition(Cell cell)
	{
	Rectangle r=cell.region.getPolygon().getBounds();
	int size=textSize;
	textFont(font, size);
	boolean fit=false;
	double diff=10000000;
	int bestSize=size;
	int bestX=-1, bestY=-1;
	
	while(!fit)
		{
		float tw=textWidth(cell.term.name);
		int y=r.y+10;
		while(y+size<r.y+r.height)
			{
			int x=r.x+5;
			while(x+tw<r.x+r.width)
				{
				Rectangle rt=new Rectangle(x-5,y-5,(int)Math.ceil(tw)+10,size+10);
				if(cell.region.getPolygon().contains(rt) &&	//here may fit
					diff>Point2D.distance((double)cell.centroid[0], (double)cell.centroid[1], (double)x+tw*0.5, (double)y+textSize*0.5))
					{
					diff=Point2D.distance((double)cell.centroid[0], (double)cell.centroid[1], (double)x+tw*0.5, (double)y+textSize*0.5);
					bestX=x; bestY=y; bestSize=size;
					fit=true;
					}
				else
					{
					x+=1;
					}
				}
			y+=1;
			}
		size--;
		
		textFont(font, size);
		if(fit==true)
			{
			cell.labelX=bestX;
			cell.labelY=bestY;
			cell.labelSize=bestSize;
			return;
			}
		else if(size<=8)
			{
			return;
			}
		}
	}
/**
 * Displaces every region in the tessellation by (x,y)
 * @param cell
 * @param x
 * @param y
 */
public void recursiveRegionTranslation(Cell cell, int x, int y)
	{
	cell.translate(x,y);
	
	//2) Continue with recursion
	if(cell.subcells!=null && cell.subcells.length>=1)//TODO: not sure why but if I do >0 the sections move...
		{
		for(Cell cc:cell.subcells)
			recursiveRegionTranslation(cc, x,y);
		}
	}

/**
 * When the mouse id double-clicked, the web browser is opened with the corresponding KEGG pathway, colored.
 */
//by gene ids
public void mouseReleased() {
    // do something based on mouse movement
	// update the screen (run draw once)
	
	if(mouseEvent.getClickCount()==1 && hoveredCells.size()>0)
		{
		if(altDown)
			{
			altDown=false;
			
			if(expData==null)	return;
			if(minHoveredCell!=null)
				{
				int x=voronto.getLocation().x+this.getWidth();
				int y=voronto.getLocation().y;
				if(heatframe!=null)	
					{
					x=heatframe.getX();
					y=heatframe.getY();
					heatframe.dispose();
					}
			
				heatframe=new JFrame();
				heatframe.setVisible(true);
				heatframe.setResizable(false);
				heatframe.setLocation(x, y);
				if(minHoveredCell.term.geneExs==null || minHoveredCell.term.geneExs.values().iterator().next().size()!=expData.getNumConditions())
					minHoveredCell.computeExpressionProfile(expData);
				gh=new CellHeatmap(minHoveredCell, font, this);
				heatframe.setTitle(minHoveredCell.term.name);
				heatframe.setSize(gh.getWidth(), gh.getHeight()+22);
				heatframe.add(gh);
				gh.init();
				gh.requestFocusInWindow();
				}
			return;
			}
		}
		if(mouseEvent.getClickCount()==2 && hoveredCells.size()>0)						//show term info on its related webpage
			{
			if(type==KEGG)
				{
				Cursor hourglassCursor = new Cursor(Cursor.WAIT_CURSOR);
				setCursor(hourglassCursor);
				
				for(Cell c:hoveredCells)
					{
					if(c.subcells==null || c.subcells.length==0)
						{
						try{
						KEGGLocator  locator = new KEGGLocator();
						KEGGPortType serv    = locator.getKEGGPort();
						
						if(expData!=null)
							{//
							PathwayElement ps[]=serv.get_elements_by_pathway("path:"+this.expData.organismKegg+c.term.id);
							if(ps.length==0)	//generic KO pathway
								{
								ps=serv.get_elements_by_pathway("path:ko"+c.term.id);//this has KO ids, no gene ids
								c.computeKOIDExpression(this.expData);
								
								long t= System.currentTimeMillis();
								String[] element_id_list=c.term.koExs.keySet().toArray(new String[0]);
								
								ArrayList<String> fg=new ArrayList<String>();
								ArrayList<String> bg=new ArrayList<String>();
								ArrayList<Integer> els=new ArrayList<Integer>();
								
								for(PathwayElement p:ps)
									{
									//System.out.println(p.getType()+"\t"+p.getElement_id());
									if(p.getType().equals("ortholog"))
										{
										float ex=0;
										int cont=0;
										
										for(String g:p.getNames())
											{
											int pos=0;
											if((pos=Arrays.binarySearch(element_id_list, g.replace("ko:", "")))>=0)
												{
												ArrayList<Float> exs=((ArrayList<Float>)c.term.koExs.get(element_id_list[pos]));
												if(exs!=null && exs.size()==this.expData.getNumConditions())
													{
													ex+=exs.get(this.selectedCol);
													cont++;
													}
												}
											}
										if(cont>0)
											{
											ex/=cont;
											switch(COLOR_MODE)
												{
												case COLOR_EXPRESSION:
													if(ex>avgExp)		
														{
														int h=(int)Math.round(255-((ex-avgExp)/(maxExp-avgExp)*255));
														bg.add("#"+Color2Hex(255, h, h));
														}
													else
														{
														int h=(int)Math.round(255-(Math.abs(ex-avgExp)/Math.abs(minExp-avgExp)*255));
														bg.add("#"+Color2Hex(h,h, 255));
														}
													break;
												case COLOR_DEVIATION:
													int h=(int)Math.round(255-(Math.abs(ex-avgExp)/Math.max(Math.abs(avgExp-minExp), Math.abs(avgExp-maxExp)))*255);
													bg.add("#"+Color2Hex(h,255,h));
													break;
												}
											fg.add("#007700");
											els.add(p.getElement_id());
											}//if it has expression mapped
										}//if path element is a gene
									}//for each path element
								int [] iels=new int[els.size()];
								for(int i=0;i<els.size();i++) iels[i]=els.get(i);
								System.out.println("Time in mapping pathway elements: "+((System.currentTimeMillis()-t)/1000.0));
								if(els.size()==0)
									{
									JOptionPane.showMessageDialog(null, "KO term "+c.term.id+" has "+c.numLeaves+" mapped "+this.expData.organism+" genes, but no dedicated pathway for that species ("+this.expData.organismKegg+")", "Error", JOptionPane.ERROR_MESSAGE);
									System.err.println("No dedicated pathway");
									}
								else
									{
									String url=serv.get_html_of_colored_pathway_by_elements("path:ko"+c.term.id, iels, fg.toArray(new String[0]), bg.toArray(new String[0]));
									System.out.println(url);
									java.awt.Desktop.getDesktop().browse(java.net.URI.create(url));
									}
								}
							else	//specific species pathway
								{
								long t= System.currentTimeMillis();
								String[] element_id_list=c.term.geneExs.keySet().toArray(new String[0]);//term.geneExs misses cap dependent things (like SPBC1683.06c)
								
								ArrayList<String> fg=new ArrayList<String>();
								ArrayList<String> bg=new ArrayList<String>();
								ArrayList<Integer> els=new ArrayList<Integer>();
								//TODO: Hay genes que se encuentran en elementos visuales del pathway pero no est‡n en term.geneExs!! 
								for(PathwayElement p:ps)
									{
								//	System.out.println(p.getType()+"\t"+p.getElement_id());
									if(p.getType().equals("gene"))
										{
										float ex=0;
										int cont=0;
										
										for(String g:p.getNames())
											{
											if(g.indexOf(":")>=0)	g=g.substring(g.indexOf(":")+1).toLowerCase();
											int pos=0;
											if((pos=Arrays.binarySearch(element_id_list, g))>=0)
												{
												ArrayList<Float> exs=((ArrayList<Float>)c.term.geneExs.get(element_id_list[pos]));
												if(exs!=null && exs.size()==this.expData.getNumConditions())
													{
													ex+=exs.get(this.selectedCol);
													cont++;
													}
												}
											}
										if(cont>0)
											{
											ex/=cont;
											switch(COLOR_MODE)
												{
												case COLOR_EXPRESSION:
													if(ex>avgExp)
														{
														int h=(int)Math.round(255-((ex-avgExp)/(maxExp-avgExp)*255));
														bg.add("#"+Color2Hex(255, h, h));
														}
													else
														{
														int h=(int)Math.round(255-(Math.abs(ex-avgExp)/Math.abs(minExp-avgExp)*255));
														bg.add("#"+Color2Hex(h,h, 255));
														}
													break;
												case COLOR_DEVIATION:
													int h=(int)Math.round(255-(Math.abs(ex-avgExp)/Math.max(Math.abs(avgExp-minExp), Math.abs(avgExp-maxExp)))*255);
													bg.add("#"+Color2Hex(h,255,h));
													break;
												}
											fg.add("#007700");
											els.add(p.getElement_id());
											}//if it has expression mapped
										}//if path element is a gene
									}//for each path element
								int [] iels=new int[els.size()];
								for(int i=0;i<els.size();i++) iels[i]=els.get(i);
								System.out.println("Time in mapping pathway elements: "+((System.currentTimeMillis()-t)/1000.0));
								if(els.size()==0)
									{
									JOptionPane.showMessageDialog(null, "KO term "+c.term.id+" has "+c.numLeaves+" mapped "+this.expData.organism+" genes, but no dedicated pathway for that species ("+this.expData.organismKegg+")", "Error", JOptionPane.ERROR_MESSAGE);
									System.err.println("No dedicated pathway");
									}
								else
									{
									String url=serv.get_html_of_colored_pathway_by_elements("path:"+this.expData.organismKegg+c.term.id, iels, fg.toArray(new String[0]), bg.toArray(new String[0]));
									System.out.println(url);
									java.awt.Desktop.getDesktop().browse(java.net.URI.create(url));
									}
								}
							}
						else	//just browse to the pathway, no map
							{
							String url="http://www.kegg.jp/kegg-bin/show_pathway?map"+c.term.id;
							System.out.println(url);
							java.awt.Desktop.getDesktop().browse(java.net.URI.create(url));
							}
							
					
						}catch(Exception e){e.printStackTrace();}
						}
					}//for each hovered cell (usually one)
				redraw();
				
				Cursor normalCursor = new Cursor(Cursor.DEFAULT_CURSOR);
				setCursor(normalCursor);
				
				}//if KEGG
			else if(type==GO)
				{
				for(Cell c:hoveredCells)
					{
					if(c.level==minHoveredLevel)
						{
						try{
							String url="http://amigo.geneontology.org/cgi-bin/amigo/term_details?term="+c.term.id;
							System.out.println("Browsing to "+url);
							java.awt.Desktop.getDesktop().browse(java.net.URI.create(url));
							}catch(Exception e){e.printStackTrace();}
						}
	
					}
				}
			else if(type==REACTOME)
				{
				
				for(Cell c:hoveredCells)
					{
					if(c.level==minHoveredLevel)
						{
						try{
							String url="http://www.reactome.org/cgi-bin/eventbrowser?DB=gk_current&ID="+c.term.id;
							System.out.println("Browsing to "+url);
							java.awt.Desktop.getDesktop().browse(java.net.URI.create(url));
							}catch(Exception e){e.printStackTrace();}
						}
					}
				}
			return;
		}//if double click
    }
    
/**
 * Converts an RGB color to its corresponding hexadecimal value
 * @param r
 * @param g
 * @param b
 * @return
 */
public String Color2Hex(int r, int g, int b)
	{
      Color c = new Color(r,g,b);
      String s=Integer.toHexString( c.getRGB() & 0x00ffffff ).toUpperCase() ;
      while(s.length()<6)
    	  s="0"+s;
      return s;
   }

/**
 * Returns the map with term as the root
 * @param term
 * @return
 */
private TreeMap<OntologyTerm, TreeMap> getMap(TreeMap<OntologyTerm, TreeMap> m, OntologyTerm term)
	{
	if(m.containsKey(term))
		return m.get(term);
	else
		{
		for(OntologyTerm ot:m.keySet())
			{
			TreeMap<OntologyTerm, TreeMap> mret=getMap(m.get(ot), term);
			if(mret!=null)	return mret;
			}
		return null;
		}
	}

/**
 * Returns the parent in the hierarchy for a given term
 * @param term
 * @param root root of the current term
 * @return
 */
private OntologyTerm searchParent(OntologyTerm term, TreeMap<OntologyTerm, TreeMap> m, OntologyTerm root)
	{
	if(m==null)	return null;
	if(m.containsKey(term))
		return root;
	else
		{
		for(OntologyTerm ot:m.keySet())
			{
			OntologyTerm parent=searchParent(term, m.get(ot), ot);
			if(parent!=null)	return parent;
			}
		return null;
		}
	}

public void keyPressed()
	{
	switch(keyCode)
		{
		case 17:	//ctrl
			ctrlDown=true;
			break;
		case 18:	//alt
			altDown=true;
			break;
			
		}
	}

public void keyReleased()
	{
	//System.out.println("key: "+keyCode);
	//TODO: check if keycodes for arrows are the same on any platform (should be)
	switch(keyCode)
		{
		case 39://right
			if(expData!=null)
				{
				selectedCol=Math.min(expData.getNumConditions()-1, selectedCol+1);
				expression2color(expData, selectedCol);
				}
			break;
		case 37://left
			selectedCol=Math.max(0, selectedCol-1);
			if(expData!=null)	expression2color(expData, selectedCol);
			break;
		case 38://up
			levelThreshold=Math.max(1, levelThreshold-1);
			setOntologyName();
			break;
		case 40://down
			if(type==VoronoiVisualization.SLIMBP || type==VoronoiVisualization.SLIMCC)	
				levelThreshold=Math.min(v.maxLevel, levelThreshold+1);
			/*else if(type==VoronoiVisualization.KEGG)		
				levelThreshold=Math.min(v.maxLevel+1, levelThreshold+1);
			else if(type==VoronoiVisualization.GO)			levelThreshold=Math.min(v.maxLevel+1,levelThreshold+1);	
			else if(type==VoronoiVisualization.REACTOME)	levelThreshold=Math.min(v.maxLevel+1,levelThreshold+1);*/
			else											
				levelThreshold=Math.min(v.maxLevel+1,levelThreshold+1);
			
			setOntologyName();
			break;
		case 10://enter
			
			searchedCells.clear();
			
			//dig deep into hierarchy
			//0) get hovered cell
			Cell cell=null;
			if(minHoveredCell!=null)	cell=minHoveredCell;
			
			//1) perform Benrhardt tessellation under it (to max 2 depth level again?)
			if(cell==null)	{System.err.println("No hovered cell"); return;}
			if(getMap(map, cell.term).size()<=1)	
				{
				System.err.println("Leaf cell, no further voronoi map");
				JOptionPane.showMessageDialog(null, "This term has zero or one subterms, tessellation will not be computed", "Warning", JOptionPane.WARNING_MESSAGE);
				return;
				}
			
			Cursor hourglassCursor = new Cursor(Cursor.WAIT_CURSOR);
			setCursor(hourglassCursor);
			
			try{tessellate(cell.term);}catch(Exception e){e.printStackTrace();}
			
			Cursor normalCursor = new Cursor(Cursor.DEFAULT_CURSOR);
			setCursor(normalCursor);
			
			break;
			
		case 8://supr
			//return back into hierarchy
			if(root.id.equals("root") && root.name.equals("root"))	return;
			searchedCells.clear();
			OntologyTerm parent=searchParent(root, map, new OntologyTerm("root", "root"));
			v=tessellations.get(parent);
			if(v==null)	try{tessellate(parent);}catch(Exception e){e.printStackTrace();}
			
			root=parent;
		 	cells=v.getCells();
			hoveredCells.clear();
			computedLabels=false;
			setOntologyName();
			
			break;
		case 17:	//ctrl
			ctrlDown=false;
			break;
		case 18:	//alt
			altDown=false;
			break;
			
		default:
			System.out.println(keyCode);
			break;
		}
	switch(key)
		{
		case KeyEvent.VK_ESCAPE:
			break;
		case 's'://change scale
			//SCALE_MODE=(SCALE_MODE+1)%3;
			SCALE_MODE=(SCALE_MODE+1)%2;	//no ontology mode
			System.out.println("SCALE_MODE: "+SCALE_MODE);
			expression2color(this.expData);
			break;
		/*case 'c'://change color
			//if(COLOR_MODE==COLOR_EXPRESSION)	COLOR_MODE=COLOR_DEVIATION;
			//else								COLOR_MODE=COLOR_EXPRESSION;
			COLOR_MODE=(COLOR_MODE+1)%2;//%3 to tests standard deviation
			expression2color(this.expData, this.selectedCol);
			break;*/
		case 'l'://cell labels
			SHOW_LABELS=!SHOW_LABELS;
			break;
		case 'p':
			JFileChooser selecFile = new JFileChooser();
			selecFile.addChoosableFileFilter(new ImageFileFilter());
			if(expData!=null)	
				{
				selecFile.setCurrentDirectory(new File(expData.filePath));
				selecFile.setSelectedFile(new File("voronto-"+expData.organismKegg+"-"+expData.conditionNames[this.selectedCol]+".tif"));
				}
			else
				selecFile.setSelectedFile(new File("voronto.tif"));
			
			int returnval = selecFile.showSaveDialog(this);

			if(returnval==JFileChooser.APPROVE_OPTION)
				{
				saving=true;
				redraw();
				try{Thread.sleep(1000);}catch(Exception e){}
				save(selecFile.getSelectedFile().getAbsolutePath());
				saving=false;
				redraw();
				}
			break;
		case 'f':	//searcher
			if(searchFrame==null)
				{
				searchFrame=new JFrame();
				searchPApplet=new SearchFrame(this);
				searchFrame.setUndecorated(true);
				
				searchPApplet.init();
				searchFrame.add(searchPApplet);
				
				searchPApplet.requestFocusInWindow();
				}
			else if(searchFrame.isVisible())	{searchFrame.setVisible(false); return;}
			
			searchFrame.setSize(searchPApplet.width, searchPApplet.height);
			searchFrame.setLocation((int)(voronto.getLocation().x+this.width*0.5-searchPApplet.width*0.5), (int)(voronto.getLocation().y+this.height*0.5-searchPApplet.height*0.5));
			searchPApplet.searchText="";
			searchPApplet.redraw();
			searchFrame.setVisible(true);
			searchPApplet.requestFocusInWindow();
			
			drawSearch=false;
			break;
		case 'd':	//difexp searcher
			if(deSearchFrame==null)
				{
				deSearchFrame=new DifExpSearch(this);
				deSearchFrame.setUndecorated(true);
				deSearchFrame.requestFocusInWindow();
				}
			else if(deSearchFrame.isVisible())	{deSearchFrame.setVisible(false); return;}
			
			deSearchFrame.setLocation((int)(voronto.getLocation().x+this.width*0.5-deSearchFrame.getWidth()*0.5), (int)(voronto.getLocation().y+this.height*0.5-deSearchFrame.getHeight()*0.5));
			deSearchFrame.setVisible(true);
			
			break;
		case 'h':	//helper
			try{
				java.awt.Desktop.getDesktop().browse(java.net.URI.create("http://vis.usal.es/~visusal/voronto/voronto/Help.html"));
			}catch(Exception e){JOptionPane.showMessageDialog(null, "Error", "Cannot open browser. Please visit http://vis.usal.es/~visusal/voronto/voronto/Help.html for help", JOptionPane.ERROR_MESSAGE);}
			break;
		case 'b': //go back to input interface
			if(voronto!=null)
				{
				if(heatframe!=null)		heatframe.dispose();
				voronto.goBack();
				}
			else
				System.err.println("Voronto is null, no chance to go back");
			break;
		case 'm'://change medium value to mean or median
			whiteValue=(whiteValue+1)%2;
			System.out.println("Changing scale mode");
			
			expression2color(this.expData);
			redraw();
			break;
		}
	redraw();
	}

public void exit()
	{
	//System.out.println("Exiting...");
	return;//i don't want esc to close the program
	}
/**
 * Computes a new tesselleation with cell as root node
 * @param cell
 */
public void tessellate(OntologyTerm term) throws Exception
	{
	boolean translate=false;
	if(tessellations.get(term)==null)	//search for the tessellation for this term as root
		{
		v=new BernhardtTessellation(getMap(map,term), expData, width, height-100, type, this.maxDepth);
		tessellations.put(term, v);
	 	translate=true;
		}
	else
		v=tessellations.get(term);
	
	//2) substitute the visualization level (the root node)
	root=term;
	cells=v.getCells();
	
	if(type==VoronoiVisualization.KEGG)	this.levelThreshold=v.maxLevel+1;
 	if(type==VoronoiVisualization.SLIMBP || type==VoronoiVisualization.SLIMCC)	this.levelThreshold=v.maxLevel;
	
	if(translate)	for(Cell c:cells)	recursiveRegionTranslation(c, 0, START_Y);
	computedLabels=false;
	if(expData!=null)	
		{
		mapExpression(expData);
		expression2color(expData);
		}
	
	if(expData!=null)	
		{
		expression2color(expData);
		for(int i= 0;i<expData.getNumConditions();i++)
			expression2color(expData,i);
		System.out.println("Finished mapping colors");
		}
	
	setOntologyName();
	hoveredCells.clear();
	return;
	}

public void mouseMoved() {
    if(v==null)	return;
	
	Cell[] cells = v.getCells();
	
	//1) Check hovered cells
	synchronized(hoveredCells)
		{
		minHoveredCell=null;
		minHoveredLevel=-1;
		hoveredCells.clear();
		
		for(int i=0; i<cells.length; i++)
			{
			hoveredCells.addAll(recursiveSearchHovered(cells[i]));
			if(hoveredCells.size()>0)	break;
			}
		}
	
	//2) If no hovered cells are selected, check if other element is selected to show help messages
	hoveredBox=-1;
	if(hoveredCells.size()==0)
		{
		if(organismBox!=null && organismBox.contains(mouseX, mouseY))
			hoveredBox=ORGANISM;
		if(conditionBox!=null && conditionBox.contains(mouseX, mouseY))
			hoveredBox=CONDITION;
		if(ontologyBox!=null && ontologyBox.contains(mouseX, mouseY))
			hoveredBox=ONTOLOGY;
		if(cellBox!=null && cellBox.contains(mouseX, mouseY))
			hoveredBox=CELL;
		if(scaleBox!=null && scaleBox.contains(mouseX, mouseY))
			hoveredBox=SCALE;
		}
	
    // update the screen (run draw once)
    redraw();
	}


/**
 * Returns a list of cells hovered by the mouse, that is, the cell hovered an all its parents.
 * @param cell
 * @return
 */
public ArrayList<Cell> recursiveSearchHovered(Cell cell)
	{
	ArrayList<Cell> ret=new ArrayList<Cell>();

	if(cell.region!=null)
		{
		//1) Check if region is hovered
		if(cell.region.getPolygon().contains(mouseX, mouseY))
			{
			if(minHoveredCell==null)	minHoveredCell=cell;
			else if(minHoveredCell.level<cell.level && levelThreshold>=cell.level)	minHoveredCell=cell;
			if(cell.level>minHoveredLevel)	minHoveredLevel=cell.level;
			ret.add(cell);
			
			//2) Continue with recursion (only if the parent is hovered)
			if(cell.subcells!=null && cell.subcells.length>0)	
				for(Cell cc:cell.subcells)
					{
					ArrayList<Cell> toadd=recursiveSearchHovered(cc);
					if(toadd!=null && toadd.size()>0)
						{
						ret.addAll(toadd);
						break;
						}
					}
			}
		}
	return ret;
	}

/**Maps expression in md (all columns) to voronoi regions
 * 
 * @param md
 */
public void expression2color(ExpressionData md)
	{
	if(md==null)
		{
		for(Cell c:cells)	
			recursiveExpressionNormalization(c, -1);
		}
	else
		{	
		for(int i= 0;i<expData.getNumConditions();i++)
			expression2color(md,i);
		}
	return;
	}

public void mapExpression(ExpressionData md)
	{
	long t=System.currentTimeMillis();
	maxSd=new float[md.getNumConditions()];
	Cell[] cells=v.getCells();
	for(Cell c:cells)
		{
		//long t1=System.currentTimeMillis();
		//System.out.println("recursive mapping for "+c.term.name+" with "+c.term.geneIds.size()+" annot genes");
		recursiveExpressionMapping(c, md);
		//System.out.println("it took "+(System.currentTimeMillis()-t1)/1000.0);
		}
	System.out.println("Mapping expression to ontology takes "+(System.currentTimeMillis()-t)/1000.0);
	return;
	}

/**Maps expression in md (in experiment col) to voronoi regions
 * 
 * @param md
 */
public void expression2color(ExpressionData md, int col)
	{
	if(expData==null)	{System.err.println("No expression data provided for mapping"); return;}
	//Different options: 
	//a) against min and max of overall md
	//b) against min and max of overall md[col] -> we will start by this one, but not sure what's best, possibly last one
	//c) against min and max of mapped cell expressions
	switch(SCALE_MODE)
		{
		case SCALE_MATRIX:
		//	System.out.println("Changing scale to mode MATRIX");
			minExp=md.min;
			maxExp=md.max;
			avgExp=md.average;
			medianExp=md.median;
			break;
		case SCALE_CONDITION:
		//	System.out.println("Changing scale to mode CONDITION");
			minExp=md.minCols[col];
			maxExp=md.maxCols[col];
			avgExp=md.averageCols[col];
			medianExp=md.medianCols[col];
			break;
		case SCALE_ONTOLOGY:
			//System.out.println("Changing scale to mode ONTOLOGY");
			minExp=1000000000;
			maxExp=-1000000000;
			avgExp=0;
			contAvgExp=0;
			break;
		}
	//System.out.println("Mapping expression for "+md.getColumnLabel(col)+" in scale mode "+SCALE_MODE);
	
	Cell[] cells=v.getCells();
	
	if(SCALE_MODE==SCALE_ONTOLOGY)
		avgExp/=contAvgExp;
	for(Cell c:cells)	
		recursiveExpressionNormalization(c, col);
	
	return;
	}

public void recursiveExpressionMapping(Cell cell, ExpressionData md)
	{
	if(cell.expressionLevel.size()==md.getNumConditions())	
		return; //already computed (in GO, we can find inside loops!)

	cell.computeExpression(md, type);
	if(cell.subcells!=null && cell.subcells.length>0)	//for internal nodes
		{
		for(Cell cc:cell.subcells)
			{
			//System.out.println("recursive mapping for "+cc.term.name);//it gets into an infinite loop in biological regulation!!!
			recursiveExpressionMapping(cc, md);
			}
		//cell.computeExpressionFromChildren(md);
		}
	/*else			//for leaf nodes
		{
		//1) Compute expression
		//long t=System.currentTimeMillis();
		cell.computeExpression(md,type);
		//System.out.println("Compute expression for "+cell.term.name+" takes "+(System.currentTimeMillis()-t)/1000.0);
		//2) Check max deviation //TODO: working in KO, not in GO. Discuss if it benefits somehow
		//for(int i=0;i<md.getNumConditions();i++)
			//if(maxSd[i]<cell.expressionDeviation.get(i))	maxSd[i]=cell.expressionDeviation.get(i);
		}*/
		
	return;
	}

public void recursiveExpressionNormalization(Cell cell, int column)
	{
	//1) Normalize color
	if(column>=0 && !Float.isNaN(cell.expressionLevel.get(column)) )
		{
		int h=-1;
		if(whiteValue==MEAN)	//raw coloring
			{
			switch(COLOR_MODE)
				{
				case COLOR_EXPRESSION:
					if(cell.expressionLevel.get(column)>avgExp)
						{
						h=(int)Math.round(255-((cell.expressionLevel.get(column)-avgExp)/(maxExp-avgExp)*255));
						cell.color.set(column,new Color(255, h, h));
						}
					else
						{
						h=(int)Math.round(255-(Math.abs(cell.expressionLevel.get(column)-avgExp)/Math.abs(minExp-avgExp))*255);
						cell.color.set(column, new Color(h,h, 255));
						}
					break;
				case COLOR_DEVIATION:
					h=(int)Math.round(255-(Math.abs(cell.expressionLevel.get(column)-avgExp)/Math.max(Math.abs(avgExp-minExp), Math.abs(avgExp-maxExp)))*255);
					cell.color.set(column, new Color(h,255,h));
					break;
				case COLOR_INTERNAL_DEVIATION://TODO: testing
					System.out.println("Getting it for "+cell.term.name+" with length "+cell.expressionDeviation.size());
					if(cell.expressionDeviation.size()>0)	
						h=(int)(255-(cell.expressionDeviation.get(column)/maxSd[column])*255);
					else	h=255;
					cell.color.set(column, new Color(h,255,h));
					break;
				}
			}
		else	//quantile coloring --> note: color scaling 
			{
			switch(COLOR_MODE)
				{
				case COLOR_EXPRESSION:
					int q=-1;
					
					if(SCALE_MODE==SCALE_CONDITION)		q=expData.getQuantile(cell.expressionLevel.get(column), column);
					else if(SCALE_MODE==SCALE_MATRIX)	q=expData.getQuantile(cell.expressionLevel.get(column));
					if(q>=50)
						{
						h=(int)Math.round(255-((q-50.0)/50)*255);
						cell.color.set(column, new Color(255, h, h));
						}
					else
						{
						h=(int)Math.round(255-((50.0-q)/50)*255);
						cell.color.set(column, new Color(h,h, 255));
						}
					break;
				case COLOR_DEVIATION:
					System.err.println("Option not supported for quantiles");
					//h=(int)Math.round(255-(Math.abs(cell.expressionLevel.get(column)-getWhiteValue())/Math.max(Math.abs(getWhiteValue()-minExp), Math.abs(getWhiteValue()-maxExp)))*255);
					//cell.color=new Color(h,255,h);
					break;
				case COLOR_INTERNAL_DEVIATION://TODO: testing
					System.err.println("Option not supported for quantiles");
					/*System.out.println("Getting it for "+cell.term.name+" with length "+cell.expressionDeviation.size());
					if(cell.expressionDeviation.size()>0)	
						{
						System.out.println(cell.expressionDeviation.get(column)+"\t"+expData.sdCols[column]);
						System.out.println((cell.expressionDeviation.get(column)/expData.sdCols[column])+"\t"+(cell.expressionDeviation.get(column)/expData.sdCols[column])*255);
						h=(int)(255-(cell.expressionDeviation.get(column)/maxSd[column])*255);
						}
					else	h=255;
					cell.color=new Color(h,255,h);
					break;*/
				}
			}
		}
	else
		{
		if(column<0)	
			{
			cell.color=new ArrayList<Color>();
			cell.color.add(0, new Color(255,255,255));
			}
		else			cell.color.set(column, new Color(255,255,255));
		}
	//2) Continue with recursion
	//if(cell.subcells!=null && cell.subcells.length>1)	
	if(cell.subcells!=null && cell.subcells.length>0)	
		for(Cell cc:cell.subcells)
			recursiveExpressionNormalization(cc, column);
	}



private class ImageFileFilter extends FileFilter {

	 /**
     * Returns the extension of a file (the three letters after the dot)
     * @param f	File to get the extension
     * @return	extension of f
     */
	public final String getExtension(File f) {
        String ext = null;
        String s = f.getName();
        int i = s.lastIndexOf('.');

        if (i > 0 &&  i < s.length() - 1) {
            ext = s.substring(i+1).toLowerCase();
        }
        return ext;
    }
	
	/**
     * Returns the extension of a string (the three letters after the dot)
     * @param s	String to get the extension
     * @return	extension of s
     */
 /*   public final String getExtension(String s) 
    	{
        String ext = null;
        int i = s.lastIndexOf('.');

        if (i > 0 &&  i < s.length() - 1)          ext = s.substring(i+1).toLowerCase();
        return ext;
    	}*/
    
    /**
     * Decides if a file is acceptable as TRN file
     * @param	f	File to check
     * @return	true if the file extension is xml or gml, false otherwise
     */
    public boolean accept(File f) 
    	{
        if (f.isDirectory()) 							       return true;

        String extension = getExtension(f);
        if (extension != null) 
        	{
        	if (extension.equals("png"))                  return true;
        	if (extension.equals("tga"))                  return true;
        	if (extension.equals("tif"))                  return true;
        	if (extension.equals("jpg"))                  return true;
        	if (extension.equals("jpeg"))                 return true;
            else							              return false;
            }

        return false;
    	}

    //The description of this filter
    /**
     * Returns the description of TRN files
     * @return	A brief description of expected files for TRN.
     */
    public String getDescription() {
        return "Image file (.png, .tga, .tif, .jpg, .jpeg)";
    }
}



public void search(String searchText) {
	// TODO Auto-generated method stub
	searchFrame.setVisible(false);
	searchPApplet.searchText="";
	System.out.println("Searching...");
	searchedCells.clear();
	
	if(searchText==null || searchText.length()==0)	{redraw(); return;}
	
	for(Cell c:cells)
		searchedCells.addAll(recursiveSearch(c, searchText, 0));
	System.out.println("Found "+searchedCells.size()+" occurrences");
	redraw();
	return;
}

/**
 * Searches terms that are differentially expressed on conditions g1 respect to conditions g2
 * Type tells if they must be upreg on g1 respect to g2 (0) or the other way round (1)
 * @param g1
 * @param g2
 * @param type
 */
public void deSearch(int[] g1, int[] g2, int type, float threshold) {
	// TODO Auto-generated method stub
	deSearchFrame.setVisible(false);
	System.out.println("Expression Search...");
	searchedCells.clear();
	
	if(g1==null || g1.length==0 || g2==null || g2.length==0)	{redraw(); return;}
	
	for(Cell c:cells)
		searchedCells.addAll(recursiveDESearch(c, g1, g2, type, threshold, 0));
	System.out.println("Found "+searchedCells.size()+" occurrences");
	redraw();
	return;
	}

public ArrayList<Cell> recursiveDESearch(Cell cell, int[] g1, int[] g2, int type, float threshold, int level)
	{
	ArrayList<Cell> retList=new ArrayList<Cell>();
	if(cell.subcells!=null && cell.subcells.length>0)
		{
		boolean add=false;
		for(Cell c:cell.subcells)
			{
			ArrayList<Cell> sublist=recursiveDESearch(c,g1, g2, type, threshold, level+1);
			retList.addAll(sublist);
			if(sublist.size()>0)	add=true;
			}
		if((add && level<=this.levelThreshold))
			retList.add(cell);
		}
	
	double e1=0;
	for(int i:g1)			e1+=cell.expressionLevel.get(i);
	double e2=0;
	for(int i:g2)			e2+=cell.expressionLevel.get(i);
	e1/=g1.length; e2/=g2.length;
	
	if(type==0 && e1>e2+threshold)	//TODO: must add some threshold
		retList.add(cell);
	else if(type==1 && e2>e1+threshold)
		retList.add(cell);
	
	return retList;
	}

public ArrayList<Cell> recursiveSearch(Cell cell, String searchText, int level)
	{
	ArrayList<Cell> retList=new ArrayList<Cell>();
	if(cell.subcells==null || cell.subcells.length==0)
		{
		if(cell.term.name.toLowerCase().contains(searchText.toLowerCase()) || cell.term.id.toLowerCase().contains(searchText.toLowerCase()))
			retList.add(cell);
		else
			for(String id:cell.term.geneIds)
				if(id.toLowerCase().contains(searchText.toLowerCase()))	{retList.add(cell); break;}
 		}
	else
		{
		boolean add=false;
		for(Cell c:cell.subcells)
			{
			ArrayList<Cell> sublist=recursiveSearch(c,searchText, level+1);
			retList.addAll(sublist);
			if(sublist.size()>0)	add=true;
			}
		if((add && level<=this.levelThreshold) || cell.term.name.toLowerCase().contains(searchText.toLowerCase()) || cell.term.id.toLowerCase().contains(searchText.toLowerCase()))
			retList.add(cell);
		}
	return retList;
	}
}
