package es.usal.voronto.model.clustering;

import java.util.ArrayList;
import java.util.Map;

import ch.usi.inf.sape.hac.experiment.DissimilarityMeasure;
import ch.usi.inf.sape.hac.experiment.Experiment;


public class ExpressionSubset implements Experiment, DissimilarityMeasure{

	
	private Map<String,ArrayList<Float>> expressions;
	public Map<String, ArrayList<Float>> getExpressions() {
		return expressions;
	}



	public void setExpressions(Map<String, ArrayList<Float>> expressions) {
		this.expressions = expressions;
	}

	private String[] names;
	
	
	public ExpressionSubset(Map<String,ArrayList<Float>> subset, String[] n)
		{
		expressions=subset;
		names=n;
		}
	
	

	@Override
	public int getNumberOfObservations() {
		return expressions.size();
	}

	@Override
	/**
	 * Euclidean distance
	 * @param experiment
	 * @param observation1
	 * @param observation2
	 * @return
	 */
	public double computeDissimilarity(Experiment experiment, int observation1,
			int observation2) {
		try{double comp=0;
		
		ArrayList<Float> profile1=null, profile2=null;
		profile1=((ExpressionSubset)experiment).getExpressions().get(names[observation1]);
		profile2=((ExpressionSubset)experiment).getExpressions().get(names[observation2]);
		
		
 		for(int i=0;i<profile1.size();i++)
			comp+=(profile1.get(i)-profile2.get(i))*(profile1.get(i)-profile2.get(i));
 		return Math.sqrt(comp);}
		catch(Exception e)
		{
			e.printStackTrace(); return -1;
		}
 		}

}
