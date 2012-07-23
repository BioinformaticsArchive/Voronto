package es.usal.voronto.view;

import java.awt.Color;
import java.awt.Toolkit;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JOptionPane;

import ch.usi.inf.sape.hac.HierarchicalAgglomerativeClusterer;
import ch.usi.inf.sape.hac.agglomeration.AgglomerationMethod;
import ch.usi.inf.sape.hac.agglomeration.WeightedAverageLinkage;
import ch.usi.inf.sape.hac.dendrogram.Dendrogram;
import ch.usi.inf.sape.hac.dendrogram.DendrogramBuilder;
import ch.usi.inf.sape.hac.dendrogram.DendrogramNode;
import ch.usi.inf.sape.hac.dendrogram.ObservationNode;
import ch.usi.inf.sape.hac.experiment.DissimilarityMeasure;

import es.usal.voronto.model.clustering.ExpressionSubset;
import es.usal.voronto.model.voronoi.Cell;
import processing.core.PApplet;
import processing.core.PFont;

/**
 * Class to display static, simple heatmaps for the genes in a given term.
 * It is invoked from VoronoiVisualization when a term is double clicked with Alt key down
 * TODO: make a previous hierarchical clustering on rows.
 * @author rodri
 *
 */
public class CellHeatmap extends PApplet
	{
	/**
	 * 
	 */
	private static final long serialVersionUID = -4416622212260260249L;
	Cell cell=null;
	int margin=10;//margin at the sides of the Frame
	int marginRows=100;//space for gene names		//TODO: compute from gene names' max width
	int marginCols=100;//space for condition names	//TODO: compute from condition names' max width
	int numRows=-1;
	int numCols=-1;
	VoronoiVisualization vv=null;
	String hoveredGene=null;
	
	PFont font;
	int size=12;//size of each expression level TODO: check if it becomes superlarge...
	private int maxItemsHeight;
	private int basicWidth;
	private ArrayList<Integer> order;
	private static int titleHeight=22;
	
	
	public CellHeatmap(Cell c, PFont font, VoronoiVisualization v)
		{
		super();
		vv=v;
		cell=c;
		if(cell.term.geneExs==null || cell.term.geneExs.size()==0)
			{
			System.err.println("No single gene expression for this term, please choose another one");
			return;
			}
		numRows=cell.term.geneExs.size();
		
		//cell.completeWithChildren();
		numCols=cell.term.geneExs.get(cell.term.geneExs.keySet().iterator().next()).size();
		
		this.font=font;

		basicWidth=width=marginRows+margin*2+size*numCols;
		this.height=marginCols+margin*2+size*numRows;
		int screenHeight=Toolkit.getDefaultToolkit().getScreenSize().height;
		maxItemsHeight=(int)Math.floor((screenHeight-marginCols-margin*2-(double)titleHeight)/size)-2;
		
		if(height>screenHeight-titleHeight)
			{
			height=marginCols+margin*2+size*maxItemsHeight;
					
			int numItems=numRows;
			while(numItems>maxItemsHeight)
				{
				width+=basicWidth;
				numItems-=maxItemsHeight;
				}
			
			if(width>Toolkit.getDefaultToolkit().getScreenSize().width)
				{
				System.err.println("Too many genes, choose a smaller term");
				JOptionPane.showMessageDialog(null, "This term has too many genes for the heatmap visualization. Please choose a smaller term.", "Too many genes", JOptionPane.INFORMATION_MESSAGE);
				setVisible(false);
				stop();
				return;
				}
			}
		this.setSize(width, height);
		computeOrder();
		}
	
	public void setup()
		{
		smooth();
		noLoop();
	 	}

	public void retrieveOrder(ArrayList<Integer> order, DendrogramNode dn)
		{
		if(dn.getLeft().getClass().toString().contains("ObservationNode"))
			order.add(((ObservationNode)dn.getLeft()).getObservation());
		else
			retrieveOrder(order, dn.getLeft());
		if(dn.getRight().getClass().toString().contains("ObservationNode"))
			order.add(((ObservationNode)dn.getRight()).getObservation());
		else
			retrieveOrder(order, dn.getRight());
		}
	
	public void computeOrder()
		{
		String[] names=cell.term.geneExs.keySet().toArray(new String[0]);
		order=new ArrayList<Integer>();
		
		if(names.length>1)
			{
			ExpressionSubset subset= new ExpressionSubset(cell.term.geneExs, names);
			DissimilarityMeasure dissimilarityMeasure = (DissimilarityMeasure)subset;
			//AgglomerationMethod agglomerationMethod = new CentroidLinkage();
			AgglomerationMethod agglomerationMethod = new WeightedAverageLinkage();
			//AgglomerationMethod agglomerationMethod = new CompleteLinkage();
			DendrogramBuilder dendrogramBuilder = new DendrogramBuilder(subset.getNumberOfObservations());
			HierarchicalAgglomerativeClusterer clusterer = new HierarchicalAgglomerativeClusterer(subset, dissimilarityMeasure, agglomerationMethod);
			clusterer.cluster(dendrogramBuilder);
			Dendrogram dendrogram = dendrogramBuilder.getDendrogram();
			
			DendrogramNode dn=dendrogram.getRoot();
			retrieveOrder(order, dn);
			}
		else
			{
			order.add(0);
			}
		}
	
	public void mouseMoved() 
		{
		redraw();
		}
	public void mouseReleased()
		{
		if(mouseEvent.getClickCount()==2 && hoveredGene!=null)
			{
			try{
			if(vv.expData.chip.equals("entrezgene"))
				java.awt.Desktop.getDesktop().browse(java.net.URI.create("http://www.ncbi.nlm.nih.gov/gene?term="+hoveredGene));
			else
				java.awt.Desktop.getDesktop().browse(java.net.URI.create("http://www.ncbi.nlm.nih.gov/gene?term="+hoveredGene+"%20AND%20"+vv.expData.organism.replace(" ", "%20")+"%5BOrganism%5D"));
			}catch(IOException e){System.out.println("Error: cannot show webpage: "+e.getMessage()); e.printStackTrace();}
			}
			
		}
	
	
	public void draw()
		{
		if(order==null)	return;
		fill(255,255,255);
		rect(0,0,width,height);
		noFill();
		
		textFont(font);
		
		//----------drawHeatmap
		stroke(154);
		
		String[] names=cell.term.geneExs.keySet().toArray(new String[0]);
		int cont=0;
		int xDisplacement=0;
		
		while(cont<numRows)
			{
			if(cont % maxItemsHeight==0)
				{
				if(cont>0)	xDisplacement+=basicWidth;
				
				//draw column labels
				fill(5);
				textAlign(LEFT, CENTER);
				for(int i=0;i<vv.expData.getNumConditions();i++)
					{
					if(mouseX>(margin+marginRows+xDisplacement+i*size) && mouseX<(margin+marginRows+xDisplacement+(i+1)*size))
						fill(0);
					else
						fill(154);
					
					pushMatrix();
					translate((float)(margin+marginRows+xDisplacement+(i+0.5)*size), marginCols);
					rotate((float)(1.5*PI));
					
					text(vv.expData.conditionNames[i],0, 0);
					
					popMatrix();
					}
				}
			
			textAlign(RIGHT, CENTER);
			
			//draw gene labels and expression levels
			int i=0;
			for(i=0;i<Math.min(numRows, maxItemsHeight) && i+cont<names.length;i++)//for each gene
				{
				String gene=null;
			
				gene=names[order.get(i+cont)];
				
				if(mouseY>(margin+marginCols+(i)*size) && mouseY<(margin+marginCols+(i+1)*size) && mouseX>xDisplacement && mouseX<xDisplacement+basicWidth)
					{
					fill(0);
					hoveredGene=gene;
					}
				else
					fill(154);
				text(gene.substring(gene.indexOf(":")+1).toUpperCase(),marginRows+xDisplacement, (float)(margin+marginCols+(i+0.5)*size));
				
				ArrayList<Float> exs=cell.term.geneExs.get(gene);
				for(int j=0;j<numCols;j++)
					{
					Color co=vv.getColor(exs, j);
					fill(co.getRed(), co.getGreen(), co.getBlue());
					
					float xcell=(float)(margin+marginRows+xDisplacement+j*size);
					float ycell=(float)(margin+marginCols+i*size);
					
					rect(xcell,ycell, size,size);
					}//for each col
				}//for each gene
			cont+=i;
			}
		
		//Draw hovered cell
		cont=0;
		xDisplacement=0;
		while(cont<numRows)
			{
			if(cont % maxItemsHeight==0)
				if(cont>0)	xDisplacement+=basicWidth;
		
			int i=0;
			for(i=0;i<Math.min(numRows, maxItemsHeight) && i+cont<names.length;i++)//for each gene
				{
				noFill();
				
				stroke(0);
				
				for(int j=0;j<numCols;j++)
					{
					float xcell=(float)(margin+marginRows+xDisplacement+j*size);
					float ycell=(float)(margin+marginCols+i*size);
					if(j==0 && (mouseY>ycell && mouseY<ycell+size) && mouseX<xcell)
						{
						strokeWeight(3);
						rect(xcell,ycell, size*numCols,size);
						strokeWeight(1);
						return;
						}
					if(i==0 && (mouseX>xcell && mouseX<xcell+size) && mouseY<ycell)
						{
						strokeWeight(3);
						rect(xcell,ycell, size,size*Math.min(maxItemsHeight, numRows-cont));
						strokeWeight(1);
						return;
						}
					if((mouseX>xcell && mouseX<xcell+size) && (mouseY>ycell && mouseY<ycell+size))
						{	
						strokeWeight(3);
						rect(xcell,ycell, size,size);
						strokeWeight(1);
						return;
						}
					}
				}
			cont+=i;
			}
			
		return;
		}
	}
