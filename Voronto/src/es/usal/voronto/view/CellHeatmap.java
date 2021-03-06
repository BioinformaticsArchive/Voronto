package es.usal.voronto.view;

import java.awt.Color;
import java.awt.Toolkit;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import ch.usi.inf.sape.hac.HierarchicalAgglomerativeClusterer;
import ch.usi.inf.sape.hac.agglomeration.AgglomerationMethod;
import ch.usi.inf.sape.hac.agglomeration.CentroidLinkage;
import ch.usi.inf.sape.hac.agglomeration.CompleteLinkage;
import ch.usi.inf.sape.hac.agglomeration.WeightedAverageLinkage;
import ch.usi.inf.sape.hac.dendrogram.Dendrogram;
import ch.usi.inf.sape.hac.dendrogram.DendrogramBuilder;
import ch.usi.inf.sape.hac.dendrogram.DendrogramNode;
import ch.usi.inf.sape.hac.dendrogram.ObservationNode;
import ch.usi.inf.sape.hac.experiment.DissimilarityMeasure;

import es.usal.voronto.model.clustering.ExpressionSubset;
import es.usal.voronto.model.expression.ExpressionData;
import es.usal.voronto.model.voronoi.Cell;
import es.usal.voronto.view.VoronoiVisualization.ImageFileFilter;
import es.usal.voronto.view.VoronoiVisualization.TextFileFilter;
import processing.core.PApplet;
import processing.core.PFont;
import processing.core.PGraphics;

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
	
	int nameType=1;
	private int scaleFactor=1;
	private ExpressionSubset subset;
	
	
	public CellHeatmap(Cell c, PFont font, VoronoiVisualization v)
		{
		super();
		vv=v;
		cell=c;
		scaleFactor=1;
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

	/* This method is not considering which side is better to add first in the 
	 * conext of other branches, it always goes first for left, then for right*/
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
	
	public void sortDendrogram(DendrogramNode dn)
		{
		if(dn.getClass().toString().contains("ObservationNode"))	
			return;
		//If both are leaf nodes, nothing to do
		if(dn.getLeft().getClass().toString().contains("ObservationNode") && dn.getRight().getClass().toString().contains("ObservationNode"))
			return;
		else
			{
			
			//Check order of the two left branches depending on the right node
			if(!dn.getLeft().getClass().toString().contains("ObervationNode"))
				{
				if(dn.getLeft().getLeft()!=null && dn.getLeft().getRight()!=null)
					{
					ArrayList<Float> rp=getAverageProfile(dn.getRight());
					ArrayList<Float> llp=getAverageProfile(dn.getLeft().getLeft());
					ArrayList<Float> lrp=getAverageProfile(dn.getLeft().getRight());
					double dll=subset.computeDissimilarity(rp, llp);
					double dlr=subset.computeDissimilarity(rp, lrp);
					if(dll<dlr)	//switch the branches
						{
						DendrogramNode temp=dn.getLeft().getRight();
						dn.getLeft().setRight(dn.getLeft().getLeft());
						dn.getLeft().setLeft(temp);
						}
					sortDendrogram(dn.getLeft());
					if(!dn.getRight().getClass().toString().contains("ObervationNode"))
						sortDendrogram(dn.getRight());
					}
				}
			//Check order of the two right branches depending on left node
			else
				{
				ArrayList<Float> lp=getAverageProfile(dn.getLeft());
				ArrayList<Float> rlp=getAverageProfile(dn.getRight().getLeft());
				ArrayList<Float> rrp=getAverageProfile(dn.getRight().getRight());
				double drl=subset.computeDissimilarity(lp, rlp);
				double drr=subset.computeDissimilarity(lp, rrp);
				if(drr<drl)	//switch the branches
					{
					DendrogramNode temp=dn.getRight().getLeft();
					dn.getRight().setLeft(dn.getRight().getRight());
					dn.getRight().setRight(temp);
					}
				sortDendrogram(dn.getRight());
				}
			}
		}

	public ArrayList<Float> getAverageProfile(DendrogramNode dn)
		{
		String[] names=cell.term.geneExs.keySet().toArray(new String[0]);
		
		if(dn.getClass().toString().contains("ObservationNode"))
			return cell.term.geneExs.get(names[((ObservationNode)dn).getObservation()]);
		else
			{
			ArrayList<Float> lp=getAverageProfile(dn.getRight());
			ArrayList<Float> rp=getAverageProfile(dn.getLeft());
			ArrayList<Float> ap=new ArrayList<Float>();
			for(int i=0;i<lp.size();i++)
				ap.add((float)((lp.get(i)+rp.get(i))/2.0));
			return ap;
			}
		}
		
	/*
	public ArrayList<Integer> retrieveOrder(ArrayList<Integer> order, DendrogramNode dn, ExpressionSubset subset)
		{
		ArrayList<Integer> membersL=new ArrayList<Integer>();
		ArrayList<Integer> membersR=new ArrayList<Integer>();
		ArrayList<Integer> members=new ArrayList<Integer>();
		if(dn.getLeft().getClass().toString().contains("ObservationNode"))
			{
			//order.add(((ObservationNode)dn.getLeft()).getObservation());
			membersL.add(((ObservationNode)dn.getLeft()).getObservation());
			}
		else
			membersL.addAll(retrieveOrder(order, dn.getLeft(), subset));
		
		if(dn.getRight().getClass().toString().contains("ObservationNode"))
			{
			//order.add(((ObservationNode)dn.getRight()).getObservation());
			membersR.add(((ObservationNode)dn.getRight()).getObservation());
			}
		else
			membersR.addAll(retrieveOrder(order, dn.getRight(), subset));
		
		if((membersR.size()==1 && membersL.size()==1) || membersR.size()>2 || membersL.size()>2)	
			{
			members.addAll(membersL);
			members.addAll(membersR);
			}
		else
			{
			if(membersL.size()==2)
				{
				double d0=subset.computeDissimilarity(subset, membersR, membersL.get(0));
				double d1=subset.computeDissimilarity(subset, membersR, membersL.get(1));
				if(d0<d1)
					{
					members.add(membersL.get(1));
					members.add(membersL.get(0));
					}
				else
					members.addAll(membersL);
				members.addAll(membersR);
				}
			else
				{
				double d0=subset.computeDissimilarity(subset, membersL, membersR.get(0));
				double d1=subset.computeDissimilarity(subset, membersL, membersR.get(1));
				members.addAll(membersL);
				if(d1<d0)
					{
					members.add(membersR.get(1));
					members.add(membersR.get(0));
					}
				else
					members.addAll(membersR);
				}
			}
		return members;
		}
	*/
	public void computeOrder()
		{
		String[] names=cell.term.geneExs.keySet().toArray(new String[0]);
		order=new ArrayList<Integer>();
		
		if(names.length>1)
			{
			subset= new ExpressionSubset(cell.term.geneExs, names);
			DissimilarityMeasure dissimilarityMeasure = (DissimilarityMeasure)subset;
			AgglomerationMethod agglomerationMethod = new CentroidLinkage();
			//AgglomerationMethod agglomerationMethod = new WeightedAverageLinkage();
			//AgglomerationMethod agglomerationMethod = new CompleteLinkage();
			DendrogramBuilder dendrogramBuilder = new DendrogramBuilder(subset.getNumberOfObservations());
			HierarchicalAgglomerativeClusterer clusterer = new HierarchicalAgglomerativeClusterer(subset, dissimilarityMeasure, agglomerationMethod);
			clusterer.cluster(dendrogramBuilder);
			Dendrogram dendrogram = dendrogramBuilder.getDendrogram();
			
			//dendrogram.dump();
			
			DendrogramNode dn=dendrogram.getRoot();
			sortDendrogram(dn);
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
			String entrezLabel=vv.expData.getSynonym(hoveredGene, vv.expData.ENTREZ);
				
			if(vv.expData.chip.equals("entrezgene"))
				java.awt.Desktop.getDesktop().browse(java.net.URI.create("http://www.ncbi.nlm.nih.gov/gene?term="+hoveredGene));
			else if(entrezLabel!=null)
				java.awt.Desktop.getDesktop().browse(java.net.URI.create("http://www.ncbi.nlm.nih.gov/gene?term="+entrezLabel));
			else
				java.awt.Desktop.getDesktop().browse(java.net.URI.create("http://www.ncbi.nlm.nih.gov/gene?term="+hoveredGene.toUpperCase()+"%20AND%20"+vv.expData.organism.replace(" ", "%20")+"%5BOrganism%5D"));
			}catch(IOException e){System.out.println("Error: cannot show webpage: "+e.getMessage()); e.printStackTrace();}
			}
			
		}
	
	public void keyReleased()
		{
		JFileChooser selecFile=null;
		int returnval=-333;
		
		switch(key)
			{
			case 'n'://change gene ids shown on cell heatmap
				nameType=(nameType+1)%3;
				redraw();
				break;
			case 'e':	//export gene ids 
				vv.export(cell.term);
				break;
			case 'p':	//print hi-res image
				selecFile = new JFileChooser();
				selecFile.addChoosableFileFilter(vv.new ImageFileFilter());
				
				if(vv.expData!=null)	
					{
					selecFile.setCurrentDirectory(new File(vv.expData.filePath));
					selecFile.setSelectedFile(new File(cell.term.name.trim()+".png"));
					}
				else
					selecFile.setSelectedFile(new File(cell.term.name.trim()+".png"));
				
				returnval = selecFile.showSaveDialog(this);

				if(returnval==JFileChooser.APPROVE_OPTION)
					{
					try{Thread.sleep(100);}catch(Exception e){}
					
					scaleFactor=3;
					PGraphics hires = createGraphics(width*scaleFactor, height*scaleFactor, JAVA2D);
					beginRecord(hires);
					hires.scale(scaleFactor);
					draw();
					endRecord();
					hires.save(selecFile.getSelectedFile().getAbsolutePath());
					scaleFactor=1;
					redraw();
					
					}
				break;

			default:
				break;
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
				
				String geneLabel=vv.expData.getSynonym(gene.toLowerCase(), nameType);
				if(geneLabel==null)	geneLabel=gene;
				//text(gene.substring(gene.indexOf(":")+1).toUpperCase(),marginRows+xDisplacement, (float)(margin+marginCols+(i+0.5)*size));
				text(geneLabel.substring(geneLabel.indexOf(":")+1).toUpperCase(),marginRows+xDisplacement, (float)(margin+marginCols+(i+0.5)*size));
				
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
