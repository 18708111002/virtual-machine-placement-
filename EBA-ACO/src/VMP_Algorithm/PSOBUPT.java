package VMP_Algorithm;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Random;


public class PSOBUPT {

	private int fatTreePod;
	public int core;
	public int aggre;
	public int edge;
	
	private int numVM;
	private int numPM;
	
	private double [][] c;

	
	private int [] vmCPU;
	private int [] vmMEM;
	
	private int [] pmCPU;
	private int [] pmMEM;
	
	private Pod[] pod;
	
	private int[] ffPlacement;
	private double totalEnergy;
	public double ff_bw;
	private int[][] throughSwitch;
	private double vmTotalTraffic;
	
	private final static double SWITCH_IDLE_ENERGY = 147 ;
	
	private final static double PM_EDGE_CAPACITY = 1000 * 0.8;
	private final static double EDGE_AGGRE_CAPACITY = 1000 * 0.8;
	private final static double AGGRE_CORE_CAPACITY = 1000 * 0.8;
	private final static double COMMUNICATION_BUFFER =  50;
	private final static double MIN_DELAY_PM_EDGE_LINK =  COMMUNICATION_BUFFER / PM_EDGE_CAPACITY;
	private final static double MIN_DELAY_EDGE_AGGRE_LINK =  COMMUNICATION_BUFFER / EDGE_AGGRE_CAPACITY;
	private final static double MIN_DELAY_AGGRE_CORE_LINK =  COMMUNICATION_BUFFER / AGGRE_CORE_CAPACITY;

	private final static double SERVER2_CPU_MIPS = 24000;
	private final static double SERVER2_MEM = 16000;
	
	private final static double HIGH_SERVER_FULL_ENERGY = 550;//252;
	private final static double SERVER_FULL_ENERGY = 420;//252;
	private final static double SERVER_IDLE_ENERGY = SERVER_FULL_ENERGY * 0.7;
	private final static double HIGH_SERVER_IDLE_ENERGY = HIGH_SERVER_FULL_ENERGY * 0.7;
	
	private final static int TYPE_PORT_10MB = 10;
	private final static int TYPE_PORT_100MB = 100;
	private final static int TYPE_PORT_1000MB = 1000;
	
	private static double _10MB_ENERGY = 0.2;
	private static double _100MB_ENERGY = 0.4;
	private static double _1000MB_ENERGY = 1.1;
	
	ArrayList<ArrayList<Object>> particles = new ArrayList<ArrayList<Object>>();
	ArrayList<ArrayList<Object>> localBest = new ArrayList<ArrayList<Object>>();
	ArrayList<Object> globalBest = new ArrayList<Object>();
	
	int[] bestPlacement = new int[numVM];
	
	private boolean[][]vParticle;
	
	private double[] fitness;
	private double bestFitness =Double.MIN_VALUE;
	private int bestIndex;
	
	
	
	
	private int numParticle;
	private boolean[] feasible;
	
	public PSOBUPT(int vm,int pm,int pnum)
	{
		this.numVM = vm;
		this.numPM = pm;
		this.numParticle = pnum;
		
		this.generateTopology(numPM);	
		this.generateThroughSwitch();
	}
	
	public int[] getPmCPU()
	{
		return this.pmCPU;
	}	
	
	public double[][] getTraffic()
	{
		return this.c;
	}
	
	public int[] getPmMEM()
	{
		return this.pmMEM;
	}
	public int[] getVmCPU()
	{
		return this.vmCPU;
	}
	public int[] getVmMEM()
	{
		return this.vmMEM;
	}
	public double[][] getVmTraffic()
	{
		return this.c;
	}
	
	
	public int[] PSOBUPT()
	{
		int[][]pop = new int[numParticle][numVM];
		
		this.generatePop(pop, numParticle);
			
		
		
		double fitness;
		
		this.initPSOBUPT(numParticle,pop);
		
//		for(int i = 0 ; i < numParticle; i++)
//		{
//			if(this.isFeasible(this.trasnTwoDim2Placement(this.particles.get(i))))
//			{
//				this.placeDetil("feasible",this.trasnTwoDim2Placement(this.particles.get(i)));	
//			}
//		}
		
		for(int i = 0 ; i < numParticle; i++)
		{
			
			fitness = this.calculateFitness(this.particles.get(i));
			
			if(fitness > this.fitness[i])
				this.copyParticle(this.particles.get(i), this.localBest.get(i));	
			if(fitness > this.bestFitness)
			{
				this.copyParticle(this.particles.get(i), this.globalBest);
			}
		}
		
		int curGen = 0;
		int terGen = 20;
		
		while(curGen < terGen)
		{
		
			for(int i = 0; i < numParticle; i++)
			{
				this.particleVelocityUpdate(this.getOneDim(this.particles.get(i)), this.vParticle[i], this.getOneDim(this.localBest.get(i)), this.getOneDim(this.globalBest));
				this.particlePositionUpdate(this.particles.get(i), this.vParticle[i], this.localBest.get(i), this.globalBest, i);
			}
			
			curGen++;
		}
		
		this.totalEnergy = this.totalEnergy(this.trasnTwoDim2Placement(this.globalBest));
		return this.trasnTwoDim2Placement(this.globalBest);
	}
	
	private void placeDetil(String name,int [] c)
	{

		int[] useCPU = new int[numPM];
		int[] useMEM = new int[numPM];
		boolean [] pmEmpty = new boolean[numPM];
		
		for(int i = 0; i < numPM; i++)
		{
			useCPU[i] = 0;
			useMEM[i] = 0;
		}
	
		
		for(int i = 0; i < numPM; i++)
			pmEmpty[i] = true;
		
		for(int i = 0; i < numVM; i++)
		{
			pmEmpty[c[i]] = false;
			System.out.println(name + ": vm"+i + " in " + " pm" + c[i]);
			useCPU[c[i]] += this.vmCPU[i];
			useMEM[c[i]] += this.vmMEM[i];
		}
		
		for(int j = 0; j < numPM; j++)
			if(pmEmpty[j])
				System.out.println(name+ ": pm" + j +" is empty");
		
		for(int i = 0 ;i < numPM; i++)
		{
			System.out.println("pm" + i + " u :" + useCPU[i] / (1.0 * this.pmCPU[i]) + 
								"   pm" + i + " m :" + useMEM[i] / (1.0 * this.pmMEM[i]));
		}
		
	}//end of placeDetil
	private boolean[] getOneDim(ArrayList<Object> particle)
	{
		boolean[] dim = new boolean[numPM];
		
		for(int i = 0; i < numPM; i++)
			dim[i] = true;
		
		for(int i = 0; i < numPM; i++)
			if(particle.get(i) instanceof Boolean)
				dim[i] = false;
		
		return dim;
	}//END OF getOneDim
	
	private void particleVelocityUpdate(boolean[] particle,boolean [] velocity,boolean [] lbest,boolean [] gbest)
	{		
		boolean [] part1 = this.subtraction(lbest, particle);
		boolean [] part2 = this.subtraction(gbest, particle);
		
		this.addition(velocity, part1, 0.5);
		this.addition(velocity, part2, 0.3);
		
	}// end of particleVelocityUpdate
	
	
	private void particlePositionUpdate(ArrayList<Object> particle,boolean[] vParticle,ArrayList<Object> lbest,ArrayList<Object> gbest,int pNo)
	{
		double fit_i = this.calculateFitness(this.particles.get(pNo));
		double fit_gbest = this.bestFitness;
		double fit_lbest = this.fitness[pNo];
		
		double p1 = this.getP1(fit_i, fit_gbest, fit_lbest);
		double p2 = this.getP2(fit_i, fit_gbest, fit_lbest);
		double p3 = this.getP3(fit_i, fit_gbest, fit_lbest);
		
		
		 ArrayList<Object> tparticle = new ArrayList<Object>();
	
			for(int j = 0; j < numPM; j++)
			{
				tparticle.add(false);
				
				if(!vParticle[j])
				{
					if(Math.random() <= p1)
					{
						Object obj = particle.get(j);
						tparticle.add(obj);
					}	
					else if(Math.random() <= p2 && Math.random() > p1)
					{
						Object obj = lbest.get(j);
						tparticle.add(obj);
						
					}
					else if(Math.random() <= p3 && Math.random() > p2)
					{
						Object obj = gbest.get(j);
						tparticle.add(obj);
					}
				}
			}
				
		if(this.isFeasible(this.trasnTwoDim2Placement(tparticle)))
		{
//			this.placeDetil("", this.trasnTwoDim2Placement(tparticle));
			double fitness = this.calculateFitness(tparticle);
			
			if(fitness > this.fitness[pNo])
			{
				this.fitness[pNo] = fitness;
				this.copyParticle(tparticle, lbest);
			}
			
			if(fitness > this.bestFitness)
			{
				this.bestFitness = fitness;
				this.copyParticle(tparticle, gbest);
			}
			
			this.copyParticle(tparticle, particle);
		}
		else
		{
			
		}
	}
	
	private boolean[] subtraction(boolean [] o1, boolean[]o2)
	{
		boolean [] result = new boolean[numPM];
		
	
			for(int j = 0; j < numPM; j++)
			{
				if(o2[j] && o1[j] == o2[j])
					result[j] = true;
				else
					result[j] = false;
			}
		
		return result;
		
	}
	
	private boolean[] addition(boolean [] o1, boolean[]o2,double p)
	{
		boolean [] updateBit = new boolean[numPM];
		
	
			for(int j = 0; j < numPM; j++)
			{
				updateBit[j] = false;
			}
		
	
			for(int j = 0; j < numPM; j++)
			{
				if(o1[j] != o2[j])
					updateBit[j] = true;
				
				if(Math.random() > p)
					o1[j] = o2[j];
		
			}
		return updateBit;
	}
	
	private double getP1(double fit_i,double fit_gbest, double fit_lbest)
	{
		return  (1/fit_i) / ( 1/fit_i + 1/fit_lbest + 1/fit_gbest);
	}// get p1
	
	private double getP2(double fit_i,double fit_gbest, double fit_lbest)
	{
		return  (1/fit_lbest) / ( 1/fit_i + 1/fit_lbest + 1/fit_gbest);
	}// get p2
	
	private double getP3(double fit_i,double fit_gbest, double fit_lbest)
	{
		return  (1/fit_gbest) / ( 1/fit_i + 1/fit_lbest + 1/fit_gbest);
	}// get p3
	
	
	private void copyParticle(ArrayList<Object> srcParticle,ArrayList<Object> desParticle)
	{
		for(int i = 0; i < numPM; i++)
		{
			if(srcParticle.get(i) instanceof Boolean)
			{
				boolean b = (boolean) srcParticle.get(i);
				desParticle.set(i, b);
			}
			else
			{
				ArrayList<Integer> srcVM = (ArrayList<Integer>) srcParticle.get(i);
				ArrayList<Integer> desVM = new ArrayList<Integer>();
			
				for(int j = 0 ; j < srcVM.size();j++)
				{
					desVM.add(srcVM.get(j));
				}
				
				desParticle.set(i, desVM);
			}
		}
	}
	
	private int[] trasnTwoDim2Placement(ArrayList<Object> particle)
	{
		int[]placement = new int[numVM];
		
		for(int i = 0 ; i < numPM;i++)
		{
			if(!(particle.get(i) instanceof Boolean))
			{
				ArrayList<Integer> vmList = (ArrayList<Integer>) particle.get(i);
				
				for(int j = 0 ; j < vmList.size();j++)
				{
					int vm = vmList.get(j);
					placement[vm] = i;
				}
			}
		}
	
		return placement;
	}// end of trasnTwoDim2Placement
	
	private double calculateFitness(ArrayList<Object> arrayList)
	{
		double fitness = 0 ;
		
		int[]c = this.trasnTwoDim2Placement(arrayList);
		
		if(this.isFeasible(c))
		{
			fitness = this.SERVER_FULL_ENERGY / (this.totalEnergy(c) );
			return fitness;
		}
		else
		{
			double sumTraffic = 0;
			
			for(int i = 0; i < numVM; i++)
				for(int j = 0; j < numVM; j++)
				{		
					sumTraffic += 4 * this.c[i][j];
				}
			
			fitness = this.SERVER_IDLE_ENERGY / (sumTraffic + 2 * this.SERVER_FULL_ENERGY * (Math.pow(fatTreePod, 3)/4));
//			System.out.println("infesible:" +fitness);
			return Double.MIN_VALUE;
		}			
	}//end of calculateFitness
	
	
	private double sumTraffic(int[] placement)
	{
		int throughSwitch;
		double sumTraffic = 0;
		
		for(int i = 0; i < numVM; i++)
			for(int j = 0; j < numVM; j++)
			{
				throughSwitch = this.throughSwitch[placement[i]][placement[j]];
				sumTraffic += throughSwitch * this.c[i][j];
			}
		
		return sumTraffic;
			
	}// end of  bandWidth
	
	private double totalEnergy(int[] c) {
		
		boolean [] pmEmpty = new boolean[numPM];
		
		for(int i = 0; i < numPM; i++)
			pmEmpty[i] = true;
		
		for(int i = 0; i < numPM; i++)
		{
			if(c[i] != -1)
				pmEmpty[c[i]] = false;
		}
		
		this.initTopology();
		
		for(int i = 0; i < numVM; i++)
			this.updateLinkLoadByVM(i, c[i], c);
		
		double pmEnergy = this.totalPmEnergy(c);
		double switchEnergy = this.totalSwitchEnergy(c);
			
//		System.out.println("GA"  +  pmEnergy);
//		System.out.println("GA"  +  switchEnergy);
		
		return pmEnergy + switchEnergy;
	}
private boolean generatePop(int[][] pop, int popSize) {
		
		int[][] useCPU = new int[popSize][numPM];
		int[][] useMEM = new int[popSize][numPM];
		boolean [] vmPlacement = new boolean[numVM];
		boolean first = true;
		boolean isLinkOverLoad = false;
		this.feasible = new boolean[popSize];
		int choose;
		int randomSize = 0;
		
		for(int i = 0; i < popSize; i++)
			this.feasible[i] = false;
		
		for(int k = 0;  k < numVM; k++ )
		{
			vmPlacement[k] = false;	
		}
		
		
		for(int i = 0; i < popSize; i++)
			for(int j = 0; j < numVM; j++)
				pop[i][j] = -1;
		
		for(int i = 0; i < popSize; i++)
			for(int j = 0; j < numPM; j++)
			{
				useCPU[i][j] = 0;
				useMEM[i][j] = 0;
			}
		
		for(int i = popSize - randomSize; i < popSize; i++)
		{
			for(int j = 0; j < numVM; j++)
			{
				pop[i][j] =  new Random().nextInt(numPM);
			}
		}
		
		for(int i = 0; i < popSize - randomSize; i++)
		{
			this.initTopology();
			for(int j = 0; j < numVM; j++)
			{
				if(first)
					choose = new Random().nextInt(numVM);
				else
				{
					choose = this.nextPlacmented(vmPlacement, numVM);
				}
				
				for(int k = 0 ; k < numPM; k++)
				{		
					if(vmPlacement[choose])
						break;
					
					for(int t = 0; t < vmPlacement.length ; t++)
					{
						if(vmPlacement[t])
						{
							int vm = t;
							int pm1 = pop[i][t];
							
							if(this.isLinkOverLoad(k, pm1, this.c[j][vm]))
							{
//								System.out.println("ga link not enough");
								isLinkOverLoad = true;
								break;
							}
						}
						
					}
		
					if(useCPU[i][k] + vmCPU[choose] < pmCPU[k] &&
					   useMEM[i][k] + vmMEM[choose] < pmMEM[k] && !isLinkOverLoad)
					{
						useCPU[i][k] += vmCPU[choose];
						useMEM[i][k] += vmMEM[choose];
						pop[i][choose] = k;
						vmPlacement[choose] = true;
						first = false;		
						
						this.updateLinkLoadByGA(choose, k, pop[i], vmPlacement);
					}
				}					
			}
			
			if(this.isFeasible(pop[i]))
				this.feasible[i] = true;
			
			for(int t = 0; t < numVM; t++)
			{				
				vmPlacement[t] = false;	
			}
			first = true;				
		}
		
		for(int i = 0; i < popSize; i++)
			if(this.feasible[i])
			{
//				this.placeDetil("feasible", pop[i]);
				return true;
			}
		
		return  false;
	}
	
private void updateLinkLoadByGA(int vm, int pm1, int[] placement,boolean [] vmPlacement) {
	// TODO Auto-generated method stub
	
	for(int i = 0; i < numVM; i++)
	{
		if(vmPlacement[i])
		{
			double traffic = this.c[i][vm];
			
			int pm2 = placement[i];

			int pod1 = this.getPodNum(pm1);
			int pod2 = this.getPodNum(pm2);
			int edge1 = this.getEdge(pod1, pm1);
			int edge2 = this.getEdge(pod2, pm2);
			int pmIndex1 = this.getPmIndex(pod1, edge1, pm1);
			int pmIndex2 = this.getPmIndex(pod2, edge2, pm2);
			
			this.pod[pod1].used = true;
			this.pod[pod2].used = true;
			this.pod[pod1].edgeUsed[edge1] = true;
			this.pod[pod2].edgeUsed[edge2] = true;
			this.pod[pod1].pLinkE[pmIndex1][edge1] += traffic;
			this.pod[pod2].pLinkE[pmIndex2][edge2] += traffic;
			
			if(pod1 == pod2 && edge1 == edge2)
				continue;
			
			int aggre1 = this.getAggre(pod1,edge1);
			int aggre2 = this.getAggre(pod2,edge2);
			
			this.pod[pod1].aggreUsed[aggre1] = true;
			this.pod[pod2].aggreUsed[aggre2] = true;
			
			this.pod[pod1].eLinkA[edge1][aggre1] += traffic;
			this.pod[pod2].eLinkA[edge2][aggre2] += traffic;
			
			if(pod1 == pod2)
				continue ;
		
			int core1 = this.getCore(pod1, aggre1);
			int core2 = this.getCore(pod2, aggre2);
			
			this.pod[pod1].aLinkC[aggre1][core1] += traffic;
			this.pod[pod2].aLinkC[aggre2][core2] += traffic;
					
		}
	}
	
}// end of updateLinkLoadByVM
	private int nextPlacmented(boolean[] vmPlacement, int numVM2) {
		
		Random random = new Random();
		int choose = random.nextInt(numVM);
	
		while(vmPlacement[choose])
		{
			choose = ++choose % numVM;
		}
		
		return choose;
	}

private boolean isFeasible(int[] c) {
		
		int [] reqCPU = new int[numPM];
	    int [] reqMEM = new int[numPM];
	 
	    for(int j = 0; j < numPM; j++)
	    {
	    	reqCPU[j] = 0;
	    	reqMEM[j] = 0;
	    }
	    
		for(int i = 0; i < numVM; i++)
		{
			if(c[i] != -1)
			{
				reqCPU[c[i]] += this.vmCPU[i];
				reqMEM[c[i]] += this.vmMEM[i];
			}
			else
			{
				return false;
			}
		}
	
		
		ArrayList<Integer> overLoadPmList =  new ArrayList<Integer>();
		ArrayList<Integer> unoverLoadPmList =  new ArrayList<Integer>();
		
		 for(int j = 0; j < numPM; j++)
			 if(reqCPU[j] > this.pmCPU[j] || reqMEM[j] > this.pmMEM[j])
			 {		
				 overLoadPmList.add(j);
			 }
		 
		for(int i = 0; i < numPM; i++)
		{
			if(!overLoadPmList.contains(i))
				unoverLoadPmList.add(i);		
		}
			
		int count = 0;
		
		while(overLoadPmList.size() > 0)
		{
			if(count ++ > this.numPM)
				break;
			
			int overLoadPm = overLoadPmList.get(new Random().nextInt(overLoadPmList.size()));
				
			ArrayList<vmRequire> overVmList = this.getVmList(overLoadPm, c);
			
				
			int vm1 = overVmList.get(new Random().nextInt(overVmList.size())).id;
				
			for(int i = 0; i < unoverLoadPmList.size(); i++)
			{
				int unLoadPm = unoverLoadPmList.get(i);
				
				if(reqCPU[unLoadPm] + this.vmCPU[vm1] < this.pmCPU[unLoadPm] && reqMEM[unLoadPm]+ this.vmMEM[vm1] < this.pmMEM[unLoadPm])
				{
					reqCPU[unLoadPm] += this.vmCPU[vm1];
					reqMEM[unLoadPm] += this.vmMEM[vm1];
					
					reqCPU[overLoadPm] -= this.vmCPU[vm1];
					reqMEM[overLoadPm] -= this.vmMEM[vm1];
					
					c[vm1] = unLoadPm;
					
					break;
				}
			}
			
			overLoadPmList = this.getOverLoadServerList(c);
		}	
		
			
		
		 for(int j = 0; j < numPM; j++)
			 if(reqCPU[j] > this.pmCPU[j] || reqMEM[j] > this.pmMEM[j])
			 {		
				 return false;
			 }
		 
		 for(int i = 0; i < numVM; i++)
			 for(int j = 0; j < numPM; j++)
			 {
				 if(this.isLinkOverLoad(c[i], c[j], this.c[i][j]))
					 return false;
			 }
		 
		return true;
	}


private ArrayList<vmRequire> getVmList(int pm, int [] placement)
{
	ArrayList<vmRequire> placementedVmList = new ArrayList<vmRequire>();
	
	for(int i = 0; i < this.numVM; i++)
		if(placement[i] == pm)
		{
			vmRequire vmr = new vmRequire(vmCPU[i],vmMEM[i],i);
			placementedVmList.add(vmr);
		}
	
	placementedVmList.sort(new SortByCPU());
	
	return placementedVmList;
}// end of getSortedVmList

private ArrayList<Integer> getOverLoadServerList(int[] c)
{
	int [] reqCPU = new int[numPM];
    int [] reqMEM = new int[numPM];
 
    for(int j = 0; j < numPM; j++)
    {
    	reqCPU[j] = 0;
    	reqMEM[j] = 0;
    }
    
	for(int i = 0; i < numVM; i++)
	{
		if(c[i] != -1)
		{
			reqCPU[c[i]] += this.vmCPU[i];
			reqMEM[c[i]] += this.vmMEM[i];
		}
		else
		{
			
		}
	}

	
	ArrayList<Integer> overLoadPmList =  new ArrayList<Integer>();
	
	 for(int j = 0; j < numPM; j++)
		 if(reqCPU[j] > this.pmCPU[j] || reqMEM[j] > this.pmMEM[j])
		 {		
			 overLoadPmList.add(j);
		 }
	
	return overLoadPmList;
	
}// end of getOverLoadServerList
	
	private void initPSOBUPT(int numParticle, int[][] pop)
	{
		this.vParticle = new boolean[numParticle][numPM];
		this.fitness = new double[numParticle];
		
		for(int i = 0; i < numParticle; i++)
		{
			this.particles.add(new ArrayList<Object>());
			this.localBest.add(new ArrayList<Object>());
			
			
			for(int j = 0 ; j < numPM; j++)
			{
				this.vParticle[i][j] = false;
				
				this.particles.get(i).add(false);
				this.localBest.get(i).add(false);
				this.globalBest.add(false);
			}
		}
		
	
		
		for(int i = 0; i < numParticle; i++)
			for(int j = 0; j  < numVM; j++)
			{
				int pm = pop[i][j];
				if(pm != -1)
				{
				
					this.vParticle[i][pm] = true;
					
					ArrayList<Integer> vmList;
					
					
					if(this.particles.get(i).get(pm) instanceof Boolean)
					{
						this.particles.get(i).set(pm, new ArrayList<Integer>());
						this.localBest.get(i).set(pm, new ArrayList<Integer>());	
					}
					
					vmList = (ArrayList<Integer>) this.particles.get(i).get(pm);
					vmList.add(j);
					
					vmList = (ArrayList<Integer>) this.localBest.get(i).get(pm);
					vmList.add(j);	
				}
			}	
	}
	
	public void getResource(int[] pmCPU, int[] pmMEM, int[] vmCPU, int[] vmMEM, double[][] c)
	{
		this.pmCPU = pmCPU;
		this.pmMEM = pmMEM;
		this.vmCPU = vmCPU;
		this.vmMEM = vmMEM;
		this.c = c;
		
		 for (int x = 0; x < this.c.length; x++) {  
	            for (int y = 0; y < this.c[x].length; y++) {  
	            	
	                this.vmTotalTraffic += this.c[x][y] ;  
	            }  
		 }
	}
	
	private void generateTopology(int numPM)
	{
		this.fatTreePod = (int) Math.ceil(StrictMath.pow(4 * numPM,1.0/3));
		
		if(this.fatTreePod % 2!=0)
			this.fatTreePod ++;
		
		this.core  = (fatTreePod / 2) * (fatTreePod / 2);
		this.aggre = (fatTreePod * fatTreePod) / 2;
		this.edge  = (fatTreePod * fatTreePod) / 2;
		
		this.pod = new Pod[fatTreePod];
		
		for(int i = 0; i < fatTreePod; i++)  // generate the fat-tree topology.
			this.pod[i] = new Pod(fatTreePod);
		
		this.core = (fatTreePod / 2) * (fatTreePod / 2);
		this.aggre = fatTreePod * fatTreePod / 2;
		this.aggre = fatTreePod * fatTreePod / 2;
		
	}//end of generateTopology.
	
	public void initTopology()
	{
		for(Pod pod : this.pod)
		{
			for(int i = 0; i < fatTreePod / 2; i++)
			{
				pod.edgeUsed[i] = false;
				pod.aggreUsed[i] = false;
				
				
				for(int j = 0; j < fatTreePod / 2; j++)
				{
					pod.eLinkA[i][j] = 0;//this.generateBackTraffic(limit);
					pod.aLinkC[i][j] = 0;//this.generateBackTraffic(limit);
					pod.pLinkE[i][j] = 0;
				}
			}
		}
	}// end of initTopology
	
	
	public int [] FF()
	{
		this.ffPlacement = new int[numVM];
		int[] useCPU = new int[numPM];
		int[] useMEM = new int[numPM];
		boolean[] vmPlancement = new boolean[numVM];
		boolean isLinkOverLoad = false;
		
		this.initTopology();
		
		for(int i = 0; i < numVM; i++)
		{
			ffPlacement[i] = -1;
			vmPlancement[i] = false;
		}
		
		for(int i = 0; i < numPM; i++)
		{
			useCPU[i] = 0;
			useMEM[i] = 0;
		}
		
		
		for(int i = 0; i < numVM; i++)
		{	
			int server = this.getServer(i, useCPU, useMEM, ffPlacement);
			
			if(server > -1)
			{
				ffPlacement[i] = server;
				this.updateLinkLoadByFF(i, server, ffPlacement);
				this.updatePmResource(i, server, useCPU, useMEM);
			}
			else
			{
				return null;
			}

		}
		
		double pmEnergy = this.totalPmEnergy(ffPlacement);
		double switchEnergy = this.totalSwitchEnergy(ffPlacement);
		
		this.totalEnergy = pmEnergy + switchEnergy;
		java.text.DecimalFormat df = new java.text.DecimalFormat("#.00");
		this.ff_bw = Double.valueOf(df.format(this.useMeanSwitch(ffPlacement)));
		
		return ffPlacement;	
	}// end of FFD.
	
	private void updatePmResource(int vm, int pm, int[] useCPU, int[] useMEM)
	{
		useCPU[pm] += this.vmCPU[vm];
		useMEM[pm] += this.vmMEM[vm];
		
	}// end of updatePmResource
	
	private int getServer(int vm, int[] useCPU, int[] useMEM, int[] placement) // return servers can hold the VM.
	{
		
		for(int i = 0; i < numPM; i++)
		{
			if(	this.isResouceSatisfy(vm,i,useCPU,useMEM) &&
				this.isNetworkSatisfy(vm,i,placement)	)
			{
				
				return i;	
			}		
		}	
		return -1;	
	}// end of getServerList
	
	private boolean isResouceSatisfy(int vm, int pm, int[] useCPU, int[] useMEM)
	{
		if( useCPU[pm] + this.vmCPU[vm] < this.pmCPU[pm] &&
			useMEM[pm] + this.vmMEM[vm] < this.pmMEM[pm] )
			return true;
		
		return false;
	}// end of isResouceSatisfy
	private boolean isNetworkSatisfy(int vm1, int pm1, int[] placement)
	{
		for(int i = 0 ; i < vm1; i++)
		{
			int vm2 = i;
			int pm2 = placement[vm2];
			
			if(this.isLinkOverLoad(pm1, pm2, this.c[vm1][vm2]))
			{
//				System.out.println("network is not enough");
				return false;
			}
			
		}
		
		return true;
	}// end of isNetworkSatisfy
	
	public double useMeanSwitch(int[] placement)
	{
		int throughSwitch;
		double sumTraffic = 0;
		
		for(int i = 0; i < numVM; i++)
			for(int j = 0; j < numVM; j++)
			{
//				System.out.println(j + " " + placement[j]);
				throughSwitch = this.throughSwitch[placement[i]][placement[j]];
				sumTraffic += throughSwitch * this.c[i][j];
			}
		
		
//		System.out.println("ffd totalTraffic:" + vmTotalTraffic);
//		System.out.println("ffd:" + sumTraffic);
		
		return sumTraffic / this.vmTotalTraffic;
			
	}// end of  bandWidth
	
	private void generateThroughSwitch()
	{
		this.throughSwitch = new int[numPM][numPM];
		
		for(int i = 0; i < numPM; i++)
			for(int j = 0 ; j < numPM; j++)
				this.throughSwitch[i][j] = this.getThroughSwitch(i,j);
		
	}// end of generateThroughSwitch
	
	private int getThroughSwitch(int pm1,int pm2)
	{
		if(pm1 == pm2)
			return 0;
		
		int pod1 = this.getPodNum(pm1);
		int pod2 = this.getPodNum(pm2);
		
		if(pod1 == pod2)
		{
			int edge1 = this.getEdge(pod1, pm1);
			int edge2 = this.getEdge(pod2, pm2);
			
			return (edge1 == edge2) ? 1 : 3;
		}
		else
			return 5;
		
	}// end of getThroughSwitch
	
	private void updateLinkLoadByFF(int vm,int pm1,int[] placement)
	{
		
		for(int i = 0; i < vm; i++)
		{
			int vm2 = i;
			double traffic = this.c[vm][vm2];
		
			int pm2 = placement[vm2];
			int pod1 = this.getPodNum(pm1);
			int pod2 = this.getPodNum(pm2);
			int edge1 = this.getEdge(pod1, pm1);
			int edge2 = this.getEdge(pod2, pm2);
			int pmIndex1 = this.getPmIndex(pod1, edge1, pm1);
			int pmIndex2 = this.getPmIndex(pod2, edge2, pm2);
			
			this.pod[pod1].used = true;
			this.pod[pod2].used = true;
			this.pod[pod1].edgeUsed[edge1] = true;
			this.pod[pod2].edgeUsed[edge2] = true;
			this.pod[pod1].pLinkE[pmIndex1][edge1] += traffic;
			this.pod[pod2].pLinkE[pmIndex2][edge2] += traffic;
			
			if(pod1 == pod2 && edge1 == edge2)
				continue;
			
			int aggre1 = this.getAggre(pod1,edge1);
			int aggre2 = this.getAggre(pod2,edge2);
			
			this.pod[pod1].aggreUsed[aggre1] = true;
			this.pod[pod2].aggreUsed[aggre2] = true;
			
			this.pod[pod1].eLinkA[edge1][aggre1] += traffic;
			this.pod[pod2].eLinkA[edge2][aggre2] += traffic;
			
			if(pod1 == pod2)
				continue ;
		
			int core1 = this.getCore(pod1, aggre1);
			int core2 = this.getCore(pod2, aggre2);
			
			this.pod[pod1].aLinkC[aggre1][core1] += traffic;
			this.pod[pod2].aLinkC[aggre2][core2] += traffic;
					
		}
		
	}
	
	private int getPodNum(int pm)
	{
		
		
		int pod = 0;
		int step = (fatTreePod * fatTreePod) / 4;
		int s = 0, e = step - 1;
		
		while(true)
		{
			if(pm >= s && pm <= e)
				break;
			else
			{
				s += step;
				e += step;
				pod ++;
			}
		}
		
		return pod;	
	}// end of getPodNum
	private int getEdge(int podNum,int pm)
	{
		int s = (fatTreePod * fatTreePod) / 4 * podNum;
		int i = pm - s;
		int step = fatTreePod / 2;
		int edge = 0;
			
		while(i >= step)
		{
			i -= step;
			edge  ++;
		}
		
		return edge;
	}//end of getEdge
	
	private int getAggre(int podNum,int edge)
	{
		double minTraffic = Double.MAX_VALUE;
		int aggre = 0;
		boolean useNewLink = true;
		
		boolean p = false;
		
		Pod pod = this.pod[podNum];
		
		for(int i = 0; i < pod.eLinkA[edge].length;i++)
		{
			if(pod.eLinkA[edge][i] < minTraffic &&  pod.eLinkA[edge][i] > 0)
			{
				useNewLink = false;
				minTraffic = pod.eLinkA[edge][i];
				aggre = i;
			}
		}
		
		if(useNewLink)
		{
			for(int i = 0; i < pod.eLinkA[edge].length;i++)
			{
				if(pod.eLinkA[edge][i] == 0 )
				{
					aggre = i;
				}
			}
		}

		return aggre;
	}//end of getAggre
	
	private int getCore(int podNum,int aggre)
	{
		double minTraffic = Double.MAX_VALUE;
		int core = 0;
		boolean useNewLink = true;
	
		Pod pod = this.pod[podNum];
		
		for(int i = 0; i < pod.aLinkC[aggre].length;i++)
		{
			if(pod.aLinkC[aggre][i] < minTraffic && pod.aLinkC[aggre][i] > 0)
			{
				useNewLink = false;
				minTraffic = pod.eLinkA[aggre][i];
				core = i;
			}
		}
		
		if(useNewLink)
		{
			for(int i = 0; i < pod.aLinkC[aggre].length;i++)
			{
				if(pod.aLinkC[aggre][i] == 0)
				{
					core = i;
				}
			}
		}
		
		return core;
	}//end of getCore
	
	private int getPmIndex(int podNum,int edge,int pm)
	{
		return ( pm - (podNum * fatTreePod * fatTreePod / 4 + edge * fatTreePod / 2));
		
	}//end of getPmIndex
	
	private boolean isLinkOverLoad(int pm1,int pm2,double traffic)
	{
		if(pm1 == pm2)
			return false;	
		
		int pod1 = this.getPodNum(pm1);
		int pod2 = this.getPodNum(pm2);
		int edge1 = this.getEdge(pod1, pm1);
		int edge2 = this.getEdge(pod2, pm2);
		
//		System.out.println("pm1:" + pm1);
//		System.out.println("pm2:" + pm2);
//
//		System.out.println("pod1:" + pod1);
//		System.out.println("pod2:" + pod2);
		
			
//		if(this.isPELinkOverLoad(edge1, pod1, pm1, traffic) || this.isPELinkOverLoad(edge2, pod2, pm2, traffic))
//			return true;
		
		int aggre1 = this.getAggre(pod1,edge1);
		int aggre2 = this.getAggre(pod2,edge2);
		
		if(this.isEALinkOverLoad(pod1, edge1, aggre1, traffic) || this.isEALinkOverLoad(pod2, edge2, aggre2, traffic))
			return true;
		
		if(this.isACLinkOverLoad(pod1, aggre1, traffic) || this.isACLinkOverLoad(pod2, aggre2, traffic))
			return true;
		
		return false;
		
		
	}// end of isLinkOverLoad
	
	private boolean isPELinkOverLoad(int edge,int podNum,int pm,double traffic)
	{
		Pod pod = this.pod[podNum];
		
		int index = pm - (podNum * fatTreePod * fatTreePod / 4 + edge * fatTreePod / 2);
		
	
		if(pod.pLinkE[index][edge] + traffic > this.PM_EDGE_CAPACITY)
		{
//			System.out.println("pm link edge overload:" + (pod.pLinkE[index][edge] + traffic) );
			return true;
		}
			
		
		return false;	
	}// end of isPELinkOverLoad
	
	private boolean isEALinkOverLoad(int podNum,int edge,int aggre,double traffic)
	{
		
		Pod pod = this.pod[podNum];
		
		if(pod.eLinkA[edge][aggre] + traffic > this.EDGE_AGGRE_CAPACITY)
		{
//			System.out.println("edge: " +edge + " aggre: " + aggre);
//			System.out.println("edge link aggre overload:" + (pod.eLinkA[edge][aggre] + traffic) );
			return true;
		}
		
		return false;
		
	}// end of isEALinkOverLoad
	
	private boolean isACLinkOverLoad(int podNum,int aggre,double traffic)
	{
		double minTraffic = Double.MAX_VALUE;
		int core = 0;
		
		Pod pod = this.pod[podNum];
		
		for(int i = 0; i < pod.aLinkC[aggre].length;i++)
		{
			if(pod.aLinkC[aggre][i] < minTraffic)
			{
				minTraffic = pod.eLinkA[aggre][i];
				core = i;
			}
		}
		
		if(pod.eLinkA[aggre][core] + traffic > this.AGGRE_CORE_CAPACITY)
		{
//			System.out.println("aggre link core overload:" + (pod.eLinkA[aggre][core] + traffic) );
			return true;
		}
		
		return false;
		
	}// end of isACLinkOverLoad
	
	private void updateLinkLoadByVM(int vm, int pm1, int[] placement) {
		// TODO Auto-generated method stub
		
		for(int i = 0; i < vm; i++)
		{
			double traffic = this.c[i][vm];
			
			int pm2 = placement[i];
			int pod1 = this.getPodNum(pm1);
			int pod2 = this.getPodNum(pm2);
			int edge1 = this.getEdge(pod1, pm1);
			int edge2 = this.getEdge(pod2, pm2);
			int aggre1 = this.getAggre(pod1,edge1);
			int aggre2 = this.getAggre(pod2,edge2);
			int pmIndex1 = this.getPmIndex(pod1, edge1, pm1);
			int pmIndex2 = this.getPmIndex(pod2, edge2, pm2);
			int core1 = this.getCore(pod1, aggre1);
			int core2 = this.getCore(pod2, aggre2);
			
			this.pod[pod1].used = true;
			this.pod[pod2].used = true;
			this.pod[pod1].edgeUsed[edge1] = true;
			this.pod[pod2].edgeUsed[edge2] = true;
			this.pod[pod1].aggreUsed[aggre1] = true;
			this.pod[pod2].aggreUsed[aggre2] = true;
			
			this.pod[pod1].pLinkE[pmIndex1][edge1] += traffic;
			this.pod[pod2].pLinkE[pmIndex2][edge2] += traffic;
			
			this.pod[pod1].eLinkA[edge1][aggre1] += traffic;
			this.pod[pod2].eLinkA[edge2][aggre2] += traffic;
			
			this.pod[pod1].eLinkA[aggre1][core1] += traffic;
			this.pod[pod2].eLinkA[aggre2][core2] += traffic;
					
		}
		
	}// end of updateLinkLoadByVM.
	
	public double totalPmEnergy(int [] placement)
	{
		int[] pmWorkCPU = new int[this.numPM];
		boolean [] pmEmpty = new boolean[numPM];
		double total = 0,u;
		
		for(int i = 0; i < numPM; i++)
		{
			pmWorkCPU[i] = 0;
			pmEmpty[i] = true;
		}
		
		for(int i = 0; i < numVM; i++)
		{
			pmEmpty[placement[i]] = false;
		}
		
		//calculate the total energy of pms.
		for(int i = 0; i < placement.length; i++)
			pmWorkCPU[placement[i]] += this.vmCPU[i]; 
				
		for(int i = 0; i < numPM; i++)
		{
			if(!pmEmpty[i])
			{
				u = (double)pmWorkCPU[i] / (double)pmCPU[i];	
				
				if(pmCPU[i] > this.SERVER2_CPU_MIPS)
					total += this.HIGH_SERVER_IDLE_ENERGY + (1 - 0.7) * this.HIGH_SERVER_FULL_ENERGY * u;
				else
					total += this.SERVER_IDLE_ENERGY + (1 - 0.7) * this.SERVER_FULL_ENERGY * u;
			}
		}
			
		return total;
	}// end of totalPmEnergy
	
	
	public double totalSwitchEnergy(int[] placement) {
		// TODO Auto-generated method stub
		double total = 0;
		double podEnergy = 0;
		
		for(int i = 0; i < this.fatTreePod; i++)
		{
			podEnergy = this.calculateEachPod(i);
//			System.out.println("pod:" + i + " " + podEnergy);
			total += podEnergy;
		}
		
		return total;
	}//end of totalSwitchEnergy
	
	private double calculateEachPod(int podNum)
	{
		Pod pod = this.pod[podNum];
		double total = 0;
		
		total += this.calcaulateEachEA(pod);
		total += this.calculateEachAC(pod);
		total += this.calcaulateEachPE(pod);
//		total += this.calcaulateEASwitch(pod);
		
		return total;
	
	}//end of calculateEachPod
	
	private double calcaulateEachEA(Pod pod)
	{
		double total = 0;
		
		for(int i = 0; i < pod.eLinkA.length;i++)
			for(int j = 0; j < pod.eLinkA[i].length; j++)
			{
				if(pod.eLinkA[i][j] < this.TYPE_PORT_10MB )
				{
						total += this._10MB_ENERGY * 2;
				}
				else if(pod.eLinkA[i][j] > this.TYPE_PORT_10MB && 
						pod.eLinkA[i][j] < this.TYPE_PORT_100MB )
				{
					total += this._100MB_ENERGY * 2;
				}
				else if(pod.eLinkA[i][j] > this.TYPE_PORT_100MB && 
						pod.eLinkA[i][j] < this.TYPE_PORT_1000MB )
				{
					total += this._1000MB_ENERGY * 2;
				}
			}
		
		return total;
	}// end of calcaulateEachEA
	
	private double calculateEachAC(Pod pod)
	{
		double total = 0;
		
		for(int i = 0; i < pod.aLinkC.length;i++)
			for(int j = 0; j < pod.aLinkC[i].length; j++)
			{
				if(pod.aLinkC[i][j] < this.TYPE_PORT_10MB )
				{
						total += this._10MB_ENERGY * 2 * 20;
				}
				else if(pod.aLinkC[i][j] > this.TYPE_PORT_10MB && 
						pod.aLinkC[i][j] < this.TYPE_PORT_100MB )
				{
					total += this._100MB_ENERGY * 2 * 20;
				}
				else if(pod.aLinkC[i][j] > this.TYPE_PORT_100MB && 
						pod.aLinkC[i][j] < this.TYPE_PORT_1000MB )
				{
					total += this._1000MB_ENERGY * 2 * 20;
				}
			}
		
		return total;
		
	}// end of calculateEachAC
	
	
	private  double calcaulateEachPE(Pod pod)
	{
		double total = 0;
		
		for(int i = 0; i < pod.pLinkE.length;i++)
			for(int j = 0; j < pod.pLinkE[i].length; j++)
			{
				if(pod.pLinkE[i][j] < this.TYPE_PORT_10MB )
				{
						total += this._10MB_ENERGY ;
				}
				else if(pod.pLinkE[i][j] > this.TYPE_PORT_10MB && 
						pod.pLinkE[i][j] < this.TYPE_PORT_100MB )
				{
					total += this._100MB_ENERGY ;
				}
				else if(pod.pLinkE[i][j] > this.TYPE_PORT_100MB && 
						pod.pLinkE[i][j] < this.TYPE_PORT_1000MB )
				{
					total += this._1000MB_ENERGY ;
				}
			}
		
		return total;
	}// end of calcaulateEachPE
	
	private double calcaulateEASwitch(Pod pod)
	{
		double total = 0;
		
		for(int i = 0; i < pod.edgeUsed.length; i++)
		{
			if(pod.edgeUsed[i])
				total += this.SWITCH_IDLE_ENERGY;
			
//			if(pod.aggreUsed[i])
//				total += this.SWITCH_IDLE_ENERGY;
		}
		
		return total;
	}// end of calcaulateEASwitch
	
	public double getTotalEnergy()
	{
		return this.totalEnergy;
	}//end of getFfdTotalEnergy

}