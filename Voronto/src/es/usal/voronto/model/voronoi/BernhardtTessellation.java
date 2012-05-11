package es.usal.voronto.model.voronoi;

import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Map;
import java.util.TreeMap;

import es.usal.voronto.control.Voronto;
import es.usal.voronto.model.expression.ExpressionData;
import es.usal.voronto.model.ontology.OntologyTerm;
import es.usal.voronto.model.voronoi.geometry.MPolygon;
import es.usal.voronto.view.VoronoiVisualization;

public class BernhardtTessellation {

	float[][] edges;//edges of the regions (hard to map to regions!)
	int width;
	int height;
	public double scale;
	private Cell[] currentCells;//current hierarchical cell structure (for computations and drawings)
	private Cell[] cells;		//total hierarchical cell structure
	int numLeaves=0;
	
	public float unit;
	private double maxDistortion=0.01;
	private boolean anyoneNeedsReweight=true;
	private Polygon boundingPolygon;
	private int contawp;
	private String currentName="root";
	
	public int maxLevel=0;//max depth ontology level found (or maxDepth if maxLevel>maxDepth)
	public int maxDepth=3;//max depth ontology level determined 
	
	private ExpressionData expData;
	
	private boolean tessellationFailed=false;
	private boolean zeroAreas=false;
	
	private int mode=0; //0-normal execution, 1-debug mode
	private int ontology;
	/**
	 * 
	 * Generates a voronoi diagram for non hierarchical cells
	 * @param cells
	 * @param width
	 * @param height
	 */
	public BernhardtTessellation(Cell[] cells, int width, int height)
		{	
		this.cells=cells;
		this.width=width;
		this.height=height;
		boundingPolygon=getBoundingBox();
		computeVoronoi();
		}
	
	

	/**
	 * TODO: Generates a voronoi diagram from an ontology hierarchy
	 * @param map
	 * @param width
	 * @param height
	 */
	public BernhardtTessellation(TreeMap<OntologyTerm,TreeMap> m, ExpressionData md, int width, int height, int ontology, int maxDepth) throws Exception
		{
		this.expData=md;
		this.ontology=ontology;
		this.maxDepth=maxDepth;
		//1) Recursively generate the cell hierarchy
		long time=System.currentTimeMillis();
		OntologyTerm[] terms=m.keySet().toArray(new OntologyTerm[0]);//TODO: this might lead to unnefficient computing.
		//cells=new Cell[terms.length];
		ArrayList<Cell> cellsTemp=new ArrayList<Cell>();
		numLeaves=0;
	    for(int i=0;i<terms.length;i++)	//Covering leaf childs (only for major nodes)
	    	{
	    	Cell c=generateRecursiveCell(m.get(terms[i]), terms[i],1);
	    	if(c!=null)	cellsTemp.add(c);
	    	}
	    //maxLevel++;//the leaf level
	    cells=cellsTemp.toArray(new Cell[0]);
	    if(cells.length==0)	throw new Exception("No mapped genes. Check that gene ids match with ontology");
	    System.out.println("Number of cells "+numLeaves);
	    System.out.println("Time in building cells "+(System.currentTimeMillis()-time)/1000.0);
		
	    //2) Sort top terms by weight (lower terms have been sorted during the recursive generation)
	    time=System.currentTimeMillis();

    	Comparator<Cell> byWeight=new Voronto.WeightComparator();
		Arrays.sort(cells, byWeight);
		
	    System.out.println("Time in sorting cells "+(System.currentTimeMillis()-time)/1000.0);
		
		//3) Set initial polygon width and height
		this.width=width;
		this.height=height;
		
		//4) Compute places recursively
		this.boundingPolygon=getBoundingBox();
		this.currentCells=cells;
		
		long t=System.currentTimeMillis();
		recursiveComputeVoronoi(cells, 0, getBoundingBox());
		System.out.println("Time in computing voronoi "+ (System.currentTimeMillis()-t)/1000.0+" seconds");
		
		return;
		}
	
	private Polygon getBoundingBox()
		{
		Polygon p=new Polygon();
		p.addPoint(0, 0);
		p.addPoint(width, 0);
		p.addPoint(width, height);
		p.addPoint(0, height);
		return p;
		}
	
	public Cell buildCell(OntologyTerm term, int level)
		{
		Cell c=new Cell(1, term, level);
		if(mode==1)	System.out.println("Creating cell "+term.name+" at level "+level);
		switch(ontology)
			{
			case VoronoiVisualization.KEGG:
				c.numLeaves=c.getNumMatchedKOTerms(expData);
				break;
			case VoronoiVisualization.GO:
			case VoronoiVisualization.SLIMBP:
			case VoronoiVisualization.SLIMCC:
				c.numLeaves=c.getNumMatchedGOTerms(expData);
				break;
			case VoronoiVisualization.REACTOME:
				c.numLeaves=c.getNumMatchedTerms(expData);
				break;
			default:
				System.err.println("No matched ontology");
				c.numLeaves=c.getNumMatchedTerms(expData);
				break;
			}
		c.weight=c.numLeaves;	//In general, the term has a size proportional to the number of genes annotated with it
		if(mode==1)	System.out.println("\t with weight "+c.weight);
		numLeaves++;
		//if(c.weight==0)	return null;			//Exception: there are gene entities and none has this term annotated
		//else			
			return c;
		}
	
	public Cell generateRecursiveCell(Map<OntologyTerm, Map> m, OntologyTerm term, int level)
		{
		OntologyTerm[] names=null;
		if(m!=null)	names=m.keySet().toArray(new OntologyTerm[0]);
		else		return null;
		if(names==null || names.length==0 || level==maxDepth)	//no children, end recursion
			{
			if(expData!=null)
				{
				Cell c=buildCell(term, level);
				if(c.weight==0)	return null; //no children and no annotations, so not shown
				else			return c;
				}
			else 
				{
				numLeaves++;
				
				if(mode==1)	System.out.println("Creating cell "+term.name+" at level "+level);
				if(term.geneIds==null || term.geneIds.size()==0)	return new Cell(1, term, level);
				else											
					{
				//	System.out.println("Adding leaf node with "+term.geneIds.size()+" genes");
					return new Cell(term.geneIds.size(), term, level);
					}
				}
			}
		else	//has children
			{
			if(level>maxLevel)
				{
				maxLevel=Math.min(maxDepth, level);
				}
			ArrayList<Cell> cs=new ArrayList<Cell>();
			float w=0;
			//NOTE: this implementation does not support that a parent has a different number of elements that the sum of the elements on its children.
			//		Therefore, if i have a node with 1 children, the node with 11 elements and the children with 9, the parent will only reflect those of the children.
			// 		Moreover, if I have a node with 29 children, but 0 annotations in the children, and 1 annotation in the parent, it won't appear 
			for(OntologyTerm n:names)
				{
				Cell sc=generateRecursiveCell((Map<OntologyTerm,Map>)(m.get(n)), n, level+1);
				if(sc!=null)
					{
					cs.add(sc);
					w+=sc.weight;
					}
				}
			//NEW
			
			Cell c=buildCell(term, level);
			if(cs.size()>0)
				{
				c.subcells=cs.toArray(new Cell[0]);
				c.weight+=w;
				c.numLeaves+=w;
				Comparator<Cell> byWeight=new Voronto.WeightComparator();
				Arrays.sort(c.subcells, byWeight);
				}
			if(c.weight>0)
				return c;
			else
				return null;
			//
			//OLD
			/*if(cs.size()>0)
				{
				Cell c=new Cell(w, term, level);
				if(mode==1)	System.out.println("else: Creating cell "+term.name+" at level "+level);
				
				c.subcells=cs.toArray(new Cell[0]);
				
				Comparator<Cell> byWeight=new Voronto.WeightComparator();
				Arrays.sort(c.subcells, byWeight);
					
				return c;
				}
			else	return null;*/
			}
		}
	
	public void recursiveComputeVoronoi(Cell[] c, int level, Polygon totalArea)
	{
//	System.out.println("Computing voronoi for level "+level);
	if(level<maxDepth)	//we only compute voronoi tessellation till maxDepth (in the case of GO complete, with 14 levels, this is adviceable)
	//if(level<2)	//we only compute voronoi tessellation till maxDepth (in the case of GO complete, with 14 levels, this is adviceable)
		{
		currentCells=c;
	    boundingPolygon=totalArea;
		contawp=0;
		if(c!=null && c.length>1)//NEW
			computeVoronoi();
		for(Cell cc:c)
			{
			if(cc.subcells!=null && cc.subcells.length>1)	
				{
				if(mode==1)	System.out.println("Computing voronoi for subcells of "+cc.term.name+" (level "+cc.level+")");
				currentName=cc.term.name;
				recursiveComputeVoronoi(cc.subcells, level+1, cc.region.getPolygon());
				}
			else if(cc.subcells!=null && cc.subcells.length==1)//TODO: a term with only one children that has a term with only one children, etc. is not reflected (only 1st level)
				{
				cc.subcells[0].region=new MPolygon(cc.region.getPoints());
				cc.subcells[0].centroid=cc.centroid;
				cc.subcells[0].position=cc.position;
				recursiveComputeVoronoi(cc.subcells, level+1, cc.region.getPolygon());
				}
			}
		}
	}

	
	public void computeVoronoi()
		{
		if(currentCells==null || currentCells.length==0 ){System.err.println("No cells defined"); return;}
		
		Cell[] cellsAnt=new Cell[currentCells.length];
		for(int i=0;i<cellsAnt.length;i++)		cellsAnt[i]=new Cell(currentCells[i]);	//deep copy
			
		if(mode==1)		if(currentName.contains("aging"))	
			System.out.println("tal");
				
		anyoneNeedsReweight=true;
		//0) Initial placement: trying to optimize space based on weights
		if(mode==1)	System.out.println("0) Initial placement for "+currentName);
		
		//gridPlacement();
		//weightedGridPlacement();
		weightedPolygonPlacement();
		
		if(mode==1)
			for(Cell c:currentCells)
				System.out.println(c.term.name+" at "+c.position[0]+", "+c.position[1]+" with weight "+c.weight+" and numLeaves "+c.numLeaves);
		
		for(int kk=0;kk<15000;kk++)
			{
			
			//1) Compute voronoi regions from initial positions: TODO: this is the point where it should be weighted voronoi regions!
			if(mode==1)	System.out.println("1) AWPVT");
			computeRegionsAWP();
		//TODO: might return with a failed tessellation, no region built!
		
	    
		  //3) Estimate error between desired area and obtained area
	      if(kk==0)//first iteration, we must also compute and save desired areas
			  {
		      float sum=0;
		      for(Cell c:currentCells)	      	sum+=c.numLeaves;
		      for(Cell c:currentCells)
		        {
		      	c.desiredArea=new MPolygon(boundingPolygon).area()*c.numLeaves/sum;
		 		}
			  }

		  	float asum=0;
		  	for(Cell c:currentCells)
				asum+=c.region.area();
		  	if(mode==1)
			  	System.out.println("Regions sum area: "+asum +"\t vs \t Bounding polygon area: "+new MPolygon(boundingPolygon).area());
		  	
			boolean unmatchingAreas=false;
		    if(new MPolygon(boundingPolygon).area()>asum+1000 || new MPolygon(boundingPolygon).area()<asum-1000)
		    	{
		    	if(mode==1)	System.err.println(currentName+": Areas do not match");
		    	unmatchingAreas=true;
		    	}
	      	
	      if(!anyoneNeedsReweight)
			  break;

	      //2) Reweight if the area is very different from the desired area
	      if(mode==1)	System.out.println("2) Reweight");
		  float ratioError=computeAreaRatios();
		  
		  if(mode==1)	 System.out.println(tessellationFailed+"\t"+ratioError+"\t"+computeRatioError(cellsAnt));
		  float errIncrease=ratioError/computeRatioError(cellsAnt);
		  
		  boolean tooMuchIncrease=errIncrease>2;
		  if(zeroAreas)	tooMuchIncrease=errIncrease>2.5;
		  
		  if(kk>0 && (unmatchingAreas || tessellationFailed || tooMuchIncrease))
      	      {
			  if(mode==1) System.out.println(currentName+": Stop & Rollback");
	    	  for(int i=0;i<cellsAnt.length;i++)		currentCells[i]=new Cell(cellsAnt[i]);	//deep copy
	    	  tessellationFailed=false;
	    	  break;
	    	  }
	      else
	      	{
	    	cellsAnt=new Cell[currentCells.length];
			for(int i=0;i<cellsAnt.length;i++)		cellsAnt[i]=new Cell(currentCells[i]);	//deep copy
	      	}
	      
	 		
	      reweighting();

	      //3) Centroidal voronoi tesselation (modifies initial positions in order to make region shapes more regular)
	      if(mode==1)	System.out.println("3) CVT");
	      lloydMethod(0.5,0);
		}//iterative thingy
			
		//F) Latest centroid computation  
		computeCentroids();
		
		if(mode==1)	System.out.println("Finishing");
		}
	

	
	public void computeRegionsAWP()
		{
		Rectangle bp=boundingPolygon.getBounds();
		if(mode==1)		System.out.println("AWP"+contawp);
		if(mode==1)		if(contawp==130)
			System.out.println("awp");
		contawp++;
		
		zeroAreas=false;
		
		AWPTesselation v=new AWPTesselation(currentCells, boundingPolygon);
		
		MPolygon[] p=v.computeEdges();
		if(p==null)
			{
			//System.err.println(currentName+": Tessellation failed");
			tessellationFailed=true;
			return;
			}
		for(MPolygon pp:p)	if(pp.getPoints()==null || pp.getPoints().size()==0)
			{
			zeroAreas=true; //return;
			}
		
		for(int i=0;i<p.length;i++)
			{
			p[i].translate(bp.x, bp.y);
			currentCells[i].region=p[i];
			}
		return;
		}
	
	/**
	 * Method to compute a centroid voronoi tessellation from the given voronoi points
	 */

	public void lloydMethod(double d, int it)
		{
		//1) Compute their centroids, which are the new generators
		computeCentroids();
		
		//2) Finish if convergence acquired, compute recursively if not
		float error=avgDistCentroidToGenerator();
		if(error<d || it>20)		
			{
//			System.out.println("Terminado, error: "+error);	
			return;
			}
		else 					
			{
			for(Cell c:currentCells)
				{
				if(this.boundingPolygon.contains(c.centroid[0], c.centroid[1]))
					c.position=c.centroid;
				}
			computeRegionsAWP();

			lloydMethod(d, it+1);
			}
		}

/**
 * Computes area ratios as described by Bernhardt et al. (different to Balzer/Deussen)
 * It is basically a measure of how different are the desired and the obtained area
 * Return the mean error on the ratios
 */
public float computeAreaRatios()
	{
	float sum=0;
	for(Cell c:currentCells)
		{	
		if(c.region.area()>0)		
			c.areaRatio=c.desiredArea/c.region.area();
		else						//c.areaRatio=c.desiredArea;
			{
			//System.err.println("Regrowing "+c.term.name);
			c.areaRatio=5.0;//In the case some area is 0, we give it a chance to grow
			//c.areaRatio=c.desiredArea;
			}
		sum+=Math.abs(1-c.areaRatio);
			
		if(mode==1)	System.out.println(c.term.name+"\tdesired: "+c.desiredArea+", actual: "+c.region.area()+"\tratio:"+c.areaRatio+" weight:"+c.weight+" leaves:"+c.numLeaves);
		
		}
	return sum/currentCells.length;
	}

/**
 * Computes the ratio error average of cells
 * 
 * @param cells
 * @return
 */
public float computeRatioError(Cell[] cells)
	{
	float err=0;
	for(Cell c:cells)
		{
		err+=Math.abs(1-c.areaRatio);
		}
	return err/=cells.length;
	}

public void reweighting()
	{
	anyoneNeedsReweight=false;
	for(Cell c:currentCells)
		{
		if((c.areaRatio-1) > maxDistortion)
			{
			c.weight=(float)(c.weight*(1+c.convergenceFactor*(c.desiredArea-c.region.area())/c.desiredArea));
			anyoneNeedsReweight=true;
			//c.convergenceFactor-=0.1;	//NOTE: not used so far
			}
		}
	}

	public Polygon getPolygon(MPolygon mp)
		{
		Polygon p=new Polygon();
		for(int i=0;i<mp.getCoords().length;i++)
			p.addPoint((int)mp.getCoords()[i][0], (int)mp.getCoords()[i][1]);
		return p;
		}

	public Cell[] getCells()
		{
		return cells;
		}

	public void computeCentroids()
		{
		for(int i=0;i<currentCells.length;i++)
		 	{
			currentCells[i].centroid=new float[2];
			currentCells[i].centroid[0]=0;
			currentCells[i].centroid[1]=0;
			float[][] cc=currentCells[i].region.getCoords();
			for(int j=0;j<cc.length;j++)
				{
				currentCells[i].centroid[0]+=cc[j][0];
				currentCells[i].centroid[1]+=cc[j][1];
				};
			currentCells[i].centroid[0]/=cc.length;
			currentCells[i].centroid[1]/=cc.length;
			}
		}

	protected boolean isEdgeShared(int face1[], int face2[]){
		for(int i = 0; i < face1.length; i++){
			int cur = face1[i];
			int next = face1[(i + 1) % face1.length];
			for(int j = 0; j < face2.length; j++){
				int from = face2[j];
				int to = face2[(j + 1) % face2.length];
				if(cur == from && next == to || cur == to && next == from)
					return true;
			}
		}
		return false;
	}

	public float avgDistCentroidToGenerator()
		{
		float dist=0;
		if(currentCells==null || currentCells.length==0)
			{
			System.err.println(currentName+": avgDistCentroidToGenerator: no generators or length different to centroids");
			return -1;
			}
		for(int i=0;i<currentCells.length;i++)
			{
			dist+=Math.sqrt((currentCells[i].position[0]-currentCells[i].centroid[0])*(currentCells[i].position[0]-currentCells[i].centroid[0])+(currentCells[i].position[1]-currentCells[i].centroid[1])*(currentCells[i].position[1]-currentCells[i].centroid[1]));
			}
		return dist/currentCells.length;
		}
	
	public void gridPlacement()
		{
		double unitWidth=boundingPolygon.getBounds().width/Math.sqrt(currentCells.length);
		double unitHeight=boundingPolygon.getBounds().height/Math.sqrt(currentCells.length);
		int x0=boundingPolygon.getBounds().x;
		int y0=boundingPolygon.getBounds().y;
		int x=x0;
		int y=y0;
		boolean fitting=false;
		while(!fitting)
			{
			fitting=true;
			x=x0;
			y=y0;
			for(int i=0;i<currentCells.length;i++)	//For each cell
				{
				currentCells[i].position[0]=-1000000;
				currentCells[i].position[1]=-1000000;
				
				while(!boundingPolygon.contains(currentCells[i].position[0], currentCells[i].position[1]))//If it's not inside the polygon, go to next cross in the grid
					{
					if(x+unitWidth*0.5>x0+boundingPolygon.getBounds().width)	//If we reach the end of a row, go to next
						{
						y+=unitHeight;
						if(y>y0+boundingPolygon.getBounds().height-unitHeight*0.5)	//If we reach the last cross in the grid, reduce scale and retry
							{
							unitWidth=0.9*unitWidth;
							unitHeight=0.9*unitHeight;
							fitting=false; 
							break;
							}
						x=x0;
						}
					currentCells[i].position[0]=(float)(x+unitWidth*0.5);
					currentCells[i].position[1]=(float)(y+unitHeight*0.5);
					x+=unitWidth;
					}
				}
			}
		return;
		}
	
	//As grid placement, but grid is divided in as many points as the sum of the weights (not the sum of the cells), and 
	//places are assigned accordingly
	public void weightedGridPlacement()
	{
	int sumWeight=0;
	for(Cell c:currentCells)	sumWeight+=c.weight;
	
	double unitWidth=boundingPolygon.getBounds().width/Math.sqrt(sumWeight);
	double unitHeight=boundingPolygon.getBounds().height/Math.sqrt(sumWeight);
	int x0=boundingPolygon.getBounds().x;
	int y0=boundingPolygon.getBounds().y;
	int x=x0;
	int y=y0;
	boolean fitting=false;
	while(!fitting)
		{
		fitting=true;
		x=x0;
		y=y0;
		for(int i=0;i<currentCells.length;i++)	//For each cell
			{
			if(mode==1)	if(i>0)	System.out.println("cell "+(i-1)+"\t"+currentCells[i-1].weight+"\tposition "+currentCells[i-1].position[0]+", "+currentCells[i-1].position[1]);
			currentCells[i].position[0]=-1000000;
			currentCells[i].position[1]=-1000000;
			float w=currentCells[i].weight;
				
			while(w>0 || !boundingPolygon.contains(currentCells[i].position[0], currentCells[i].position[1]))//If it's not inside the polygon, go to next cross in the grid 
				{																						
				if(x+unitWidth*0.5>x0+boundingPolygon.getBounds().width)	//If we reach the end of a row, go to next
					{
					y+=unitHeight;
					if(y+unitHeight*0.5>y0+boundingPolygon.getBounds().height)	//If we reach the last cross in the grid, reduce scale and retry
						{
						if(mode==1)	System.out.println("cell "+i+" need to refactor, unitW and H"+unitWidth+", "+unitHeight);
						unitWidth=0.99*unitWidth;
						unitHeight=0.99*unitHeight;
						fitting=false; 
						break;
						}
					x=x0;
					}
				currentCells[i].position[0]=(float)(x+unitWidth*0.5);
				currentCells[i].position[1]=(float)(y+unitHeight*0.5);
				x+=unitWidth;
				w--;
				}
			 if(!boundingPolygon.contains(currentCells[i].position[0], currentCells[i].position[1]))//if the position falls outside the polygon
			 	{
				unitWidth=0.99*unitWidth;
				unitHeight=0.99*unitHeight;
				fitting=false; 
				break; 
			 	}
			}//for each cell
		}
	return;
	}
	
	public void weightedPolygonPlacement()
		{
		int sumWeight=0;
		for(Cell c:currentCells)	sumWeight+=c.weight;
		if(mode==1)
			if(currentName.contains("cellular process"))
				System.out.println("response to sitmulus");
		double unitWidth=boundingPolygon.getBounds().width/Math.sqrt(sumWeight);
		double unitHeight=boundingPolygon.getBounds().height/Math.sqrt(sumWeight);
		if(unitWidth==0 || unitHeight==0)
			{
			//System.err.println("Empty starting area");
			return;
			}
		int x0=boundingPolygon.getBounds().x;
		int y0=boundingPolygon.getBounds().y;
		boolean fitting=false;
		ArrayList<Point2D.Float> inits=new ArrayList<Point2D.Float>();
		inits.add(new Point2D.Float(x0,y0));
		
		float sumArea=0;
		
		
		while(!fitting)
			{
			inits.clear();
			inits.add(new Point2D.Float(x0,y0));
			fitting=true;
			for(Cell c: currentCells)
				sumArea+=Math.sqrt(c.weight)*unitWidth*Math.sqrt(c.weight)*unitHeight;
			if(mode==1)
				{
				System.out.println("Unit width and height: "+unitWidth+", "+unitHeight);
				System.out.println("Bounds: "+boundingPolygon.getBounds().width+", "+boundingPolygon.getBounds().height);
				System.out.println("Total area is "+boundingPolygon.getBounds().width*boundingPolygon.getBounds().height);
				System.out.println("Total polygonal area is "+new MPolygon(boundingPolygon).area());
				System.out.println("Total subcells area is "+sumArea);
				}
			sumArea=0;
		
			for(int i=0;i<currentCells.length;i++)	//For each cell
				{
				currentCells[i].position[0]=-1000000;
				currentCells[i].position[1]=-1000000;
				float w=currentCells[i].weight;
				Point2D.Float center=new Point2D.Float(-1000000,-1000000);
				if(mode==1)					System.out.println("Placing "+i+" with weight "+w);
				
				//int cont=0;
				//while(cont < inits.size())//keep trying
				while(inits.size()>0)
					{
					//Point2D.Float init=inits.get(cont);//get one possible init
					Point2D.Float init=inits.remove(0);//get one possible init
					
					Rectangle2D.Float r=new Rectangle2D.Float(init.x,init.y,(int)(Math.sqrt(w)*unitWidth),(int)(Math.sqrt(w)*unitHeight));
					center=new Point2D.Float(init.x+(int)(Math.sqrt(w)*unitWidth*.5), init.y+(int)(Math.sqrt(w)*unitHeight*.5));
					
					//if fitting, because it's occupied, if not, because it's out of boundingPolygon, we add the square and inits
					//inits.remove(cont);
					//add up to two new inits
					Point2D.Float np=new Point2D.Float(init.x+r.width, init.y);
					if(init.x+r.width<x0+boundingPolygon.getBounds().width && !inits.contains(np))
						inits.add(np);
						
					np=new Point2D.Float(init.x, init.y+r.height);
					if(init.y+r.height<y0+boundingPolygon.getBounds().height)
						inits.add(np);
						
					if(boundingPolygon.contains(center))	//this one is good, put position	
						{
						currentCells[i].position[0]=center.x;
						currentCells[i].position[1]=center.y;
						
						sumArea+=r.width*r.height;
						//TODO: check that the added polygon did not swallow some inits
						for(int k=0;k<inits.size();k++)
							{
							if(r.contains(inits.get(k)))	
								{
								if(mode==1) System.out.println("point "+k+" swallowed, removing");
								inits.remove(k);
								}
							}
						
						break;
						}
					
					
					//cont++;
					}
				
				//if(cont>inits.size())	//we ran out of positions, need to rescale
				if(inits.size()==0)
					{
					if(mode==1) System.out.println("Rescaling, not fitted at "+i+" with "+unitWidth+", "+unitHeight);
					unitWidth*=0.99;
					unitHeight*=0.99;
					
					fitting=false;
					i=currentCells.length;
					break;
					}
				}
			}
		return;
		}
	/*
	public void weightedGridPlacement()
	{
	int sumWeight=0;
	for(Cell c:currentCells)	sumWeight+=c.weight;
	
	double unitWidth=boundingPolygon.getBounds().width/Math.sqrt(sumWeight);
	double unitHeight=boundingPolygon.getBounds().height/Math.sqrt(sumWeight);
	int x0=boundingPolygon.getBounds().x;
	int y0=boundingPolygon.getBounds().y;
	int x=x0;
	int y=y0;
	boolean fitting=false;
	while(!fitting)
		{
		fitting=true;
		x=x0;
		y=y0;
		for(int i=0;i<currentCells.length;i++)	//For each cell
			{
			if(mode==1)	if(i>0)	System.out.println("cell "+(i-1)+"\t"+currentCells[i-1].weight+"\tposition "+currentCells[i-1].position[0]+", "+currentCells[i-1].position[1]);
			currentCells[i].position[0]=-1000000;
			currentCells[i].position[1]=-1000000;
			float w=currentCells[i].weight;
				
			while(w>0 || !boundingPolygon.contains(currentCells[i].position[0], currentCells[i].position[1]))//If it's not inside the polygon, go to next cross in the grid 
					{																						//We added w>0 to increase also as many times as its weight
					if(x+unitWidth*0.5>x0+boundingPolygon.getBounds().width)	//If we reach the end of a row, go to next
						{
						y+=unitHeight;
						if(y+unitHeight*0.5>y0+boundingPolygon.getBounds().height)	//If we reach the last cross in the grid, reduce scale and retry
							{
							//System.out.println("cell "+i+" need to refactor, unitW and H"+unitWidth+", "+unitHeight);
							unitWidth=0.99*unitWidth;
							unitHeight=0.99*unitHeight;
							fitting=false; 
							break;
							}
						x=x0;
						}
					currentCells[i].position[0]=(float)(x+unitWidth*0.5);
					currentCells[i].position[1]=(float)(y+unitHeight*0.5);
					x+=unitWidth;
					//w--;
					}
				}
			}
	return;
	}
	*/
}