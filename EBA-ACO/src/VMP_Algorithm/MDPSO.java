package VMP_Algorithm;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;


public class MDPSO {

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
	
	private int[] psoysPlacement;
	private double psoys_energy;
	public double psoys_bw;
	
	private int[][] throughSwitch;
	private double vmTotalTraffic;
	
	private ArrayList< ArrayList<Integer>> pmContains = new ArrayList<ArrayList<Integer>>();
	
	private final static double SERVER2_CPU_MIPS = 24000;
	private final static double SERVER2_MEM = 16000;
	
	private final static double HIGH_SERVER_FULL_ENERGY = 550;//252;
	private final static double SERVER_FULL_ENERGY = 420;//252;
	private final static double SERVER_IDLE_ENERGY = SERVER_FULL_ENERGY * 0.7;
	private final static double HIGH_SERVER_IDLE_ENERGY = HIGH_SERVER_FULL_ENERGY * 0.7;
	
	
	private final static double SWITCH_IDLE_ENERGY = 147 ;
	
	private final static double PM_EDGE_CAPACITY = 1000 * 0.8;
	private final static double EDGE_AGGRE_CAPACITY = 1000 * 0.8;
	private final static double AGGRE_CORE_CAPACITY = 1000 * 0.8;
	private final static double COMMUNICATION_BUFFER =  50;
	private final static double MIN_DELAY_PM_EDGE_LINK =  COMMUNICATION_BUFFER / PM_EDGE_CAPACITY;
	private final static double MIN_DELAY_EDGE_AGGRE_LINK =  COMMUNICATION_BUFFER / EDGE_AGGRE_CAPACITY;
	private final static double MIN_DELAY_AGGRE_CORE_LINK =  COMMUNICATION_BUFFER / AGGRE_CORE_CAPACITY;

	
	private final static int TYPE_PORT_10MB = 10;
	private final static int TYPE_PORT_100MB = 100;
	private final static int TYPE_PORT_1000MB = 1000;
	
	private static double _10MB_ENERGY = 0.2;
	private static double _100MB_ENERGY = 0.4;
	private static double _1000MB_ENERGY = 1.1;
	
	private int numParticle;
	boolean[][][] x_position ;
	boolean[][][] v_velocity ;
	boolean[][][] lBest;
	boolean[][] gBest;
	private boolean[] feasible;
	private double[] fitness;
	private int worstIndex;
	private int[] bestIndividual;
	private double bestFitness = Double.MAX_VALUE;
	
	private int [][] localBest; 
	private int [] globalBest;
	

	
	
	public MDPSO(int vm,int pm,int p)
	{
		this.numVM = vm;
		this.numPM = pm;
		this.numParticle = p;
		
		this.localBest = new int[numParticle][numVM];
		this.globalBest = new int [numVM];
		
		this.generateTopology(numPM);	
		this.generateThroughSwitch();
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
	
	
	public int [] MDPSO()
	{
		FFA evmp = new FFA(numVM,numPM);
		evmp.getResource(getPmCPU(), getPmMEM(), getVmCPU(), getVmMEM(), getTraffic());
		
		int[][]particles = new int[numParticle][numVM];
		
		int[][]velocity_1 = new int[numParticle][numVM];
		int[][]velocity_2 = new int[numParticle][numVM];
		
		int []t1 = new int[numVM];
		int []t2 = new int[numVM];
		
		this.initMDPSO(evmp.FFA());
		double [] fitness = new double [numParticle];
		
		for(int i = 0; i < numParticle; i++)
		{
			fitness[i] = this.particleFitness(particles[i]);
			System.arraycopy(particles[i], 0,this.localBest[i], 0, numVM);
			
			if(fitness[i] < this.bestFitness)
			{
				this.bestFitness = fitness[i];
				System.arraycopy(this.localBest[i], 0,this.globalBest, 0, numVM);
			}
		}
		
		int curGen = 0;
		int terGen = 150;
		
		while(curGen < terGen)
		{
			for(int i = 0 ;i < this.numParticle; i++)
			{
				this.updateVelocity(velocity_1[i],velocity_2[i],particles[i],this.localBest[i]);
			}
			
			curGen ++;
		}
		
		return this.globalBest;
		
	}// end of MDPSO
	
	private void updateVelocity(int [] velocity_1,int [] velocity_2,int [] particle,int []local)
	{
		int []t1 = new int[numVM];
		int []t2 = new int[numVM];
		
		double w = Math.random();
		this.multiplicationOperator(w, velocity_1, velocity_2);
		
		this.subtractionOperator(local, particle, t1, t2);
		this.multiplicationOperator(w, t1, t2);
		
		this.subtractionOperator(globalBest, particle, t1, t2);
		this.multiplicationOperator(w, t1, t2);
		
		System.arraycopy(t1, 0,velocity_1, 0, numVM);
		System.arraycopy(t2, 0,velocity_2, 0, numVM);
		
		
	}// end of updateVelocity
	
	private void transitionOperator(int [] particle,int [] velocity_1,int [] velocity_2)
	{
		for(int i = 0; i < this.numVM; i++)
		{
			if(particle[i] == velocity_1[i])
				particle[i] = velocity_2[i];
		}
		
	}// end of transitionOperator
	
	private void subtractionOperator(int [] particle1,int [] particle2,int []t1,int []t2)
	{
		for(int i = 0 ; i < this.numVM; i++)
		{
			t1[i] = particle1[i];
			t2[i] = particle2[i];
		}
		
	}// end of subtractionOperator
	
	private void additionOperator(int []t1,int []t2,int [] velocity_1_i,int [] velocity_2_i,int [] velocity_1_j,int [] velocity_2_j)
	{
		for(int i = 0 ; i < this.numVM; i++)
		{
			if(velocity_2_i[i] == velocity_1_j[i])
			{
				t1[i] = velocity_1_i[i];
				t2[i] = velocity_2_j[i];
			}
			else
			{
				t1[i] = velocity_1_i[i];
				t2[i] = velocity_2_i[i];
			}
		}
	}// end of additionOperator
	
	private void multiplicationOperator(double c ,int [] velocity_1,int [] velocity_2 )
	{
		for(int i = 0; i < this.numVM; i++)
		{
			if(Math.random() <= c)
			{
				velocity_2[i] = velocity_1[i];	
			}
		}
		
	}// end of multiplicationOperator
	
	
	private double particleFitness(int [] placement)
	{
		double tFitness = 0;
		double fitness = 0;
		
		for(int i = 0; i < numVM; i++)
		{
			int pm = placement[i];
			
			tFitness = this.pmEnergy(this.getPmUtilization(pm,this.pmContains.get(pm)),pm);
			
			this.pmContains.get(pm).remove(i);
			
			tFitness -= this.pmEnergy(this.getPmUtilization(pm,this.pmContains.get(pm)),pm);
		
			this.pmContains.get(pm).add(i);
			
			fitness += tFitness;
		}
		
		return tFitness;
		
	}// end of particleFitness
	
	private double pmEnergy(double u,int pm)
	{
		double total = 0;
				
		if(pmCPU[pm] > this.SERVER2_CPU_MIPS)
			total += this.HIGH_SERVER_IDLE_ENERGY + (1 - 0.7) * this.HIGH_SERVER_FULL_ENERGY * u;
		else
			total += this.SERVER_IDLE_ENERGY + (1 - 0.7) * this.SERVER_FULL_ENERGY * u;
		
		return total;
	}// end of pmEnergy
	
	private double getPmUtilization(int pm, ArrayList<Integer> vmList)
	{
		double u = 0 ;
		double useCPU = 0;
//		double useMEM = 0;
		
		for(int i = 0 ; i < vmList.size(); i++)
		{
			int vm = vmList.get(i);
			useCPU += this.vmCPU[vm];
//			useMEM += this.vmMEM[vm];
			
		}//end of for.
		
		u = useCPU / this.pmCPU[pm];
		
		return u;
		
	}// end of getPmUtilization
	private void initMDPSO(int [] placement)
	{
		for(int i = 0; i < numPM; i++)
		{
			this.pmContains.add(new ArrayList<Integer>());
		}
		
		
		Random random = new Random();
		
		for(int i = 0; i < numVM; i++)
		{
			int pm = random.nextInt(numPM);
			placement[i] = pm;
			
			this.pmContains.get(pm).add(i);
		}
		
	}// end of initParameter
	public int [] PSOYS()
	{
		int[][]pop = new int[numParticle][numVM];
		this.lBest = new boolean [numParticle][numVM][numPM];
		this.gBest = new boolean [numVM][numPM];
		this.x_position = new boolean[numParticle][numVM][numPM];
		this.v_velocity = new boolean[numParticle][numVM][numPM];
		int curGen = 0;
		int terGen = 20;
		
		while(curGen < terGen)
		{
			for(int i = 0; i < this.numParticle; i++)
			{
				this.generatePop(pop, numParticle);	
				this.initParticlesPop(pop);
				this.particleVelocityUpdate(this.x_position[i], this.v_velocity[i], this.lBest[i], this.gBest);
				this.particlePositionUpdate(this.x_position[i], this.v_velocity[i], this.lBest[i], this.gBest,i);
			}
			
			curGen ++;
		}
		
		return this.matrixToArry(this.gBest);
		
	}// end of PsoYsEVMP.
	
	private boolean[][] subtraction(boolean [][] o1, boolean[][]o2)
	{
		boolean [][] result = new boolean[numVM][numPM];
		
		for(int i = 0; i < numVM; i++)
			for(int j = 0; j < numPM; j++)
			{
				if(o2[i][j] && o1[i][j] == o2[i][j])
					result[i][j] = true;
				else
					result[i][j] = false;
			}
		
		return result;
		
	}
	
	private boolean[][] addition(boolean [][] o1, boolean[][]o2,double p)
	{
		boolean [][] updateBit = new boolean[numVM][numPM];
		
		for(int i = 0; i < numVM; i++)
			for(int j = 0; j < numPM; j++)
			{
				updateBit[i][j] = false;
			}
		
		for(int i = 0; i < numVM; i++)
			for(int j = 0; j < numPM; j++)
			{
				if(o1[i][j] != o2[i][j])
					updateBit[i][j] = true;
				
				if(Math.random() > p)
					o1[i][j] = o2[i][j];
		
			}
		return updateBit;
	}
	private void multiplication()
	{
		
	}
	private void particleVelocityUpdate(boolean [][] particle,boolean [][] velocity,boolean [][] lbest,boolean [][] gbest)
	{		
		boolean [][] part1 = this.subtraction(lbest, particle);
		boolean [][] part2 = this.subtraction(gbest, particle);
		
		this.addition(velocity, part1, 0.5);
		this.addition(velocity, part2, 0.3);
		
	}// end of particleVelocityUpdate
	
	private void particlePositionUpdate(boolean [][] particle,boolean [][] velocity,boolean [][] lbest,boolean [][] gbest,int index)
	{
		double fit_i = this.fitness(this.matrixToArry(particle));
		double fit_gbest = this.fitness(this.matrixToArry(gbest));
		double fit_lbest = this.fitness(this.matrixToArry(lbest));
		
		double p1 = this.getP1(fit_i, fit_gbest, fit_lbest);
		double p2 = this.getP2(fit_i, fit_gbest, fit_lbest);
		double p3 = this.getP3(fit_i, fit_gbest, fit_lbest);
		
		boolean[][] tparticle = new boolean[numVM][numPM];
		
		this.copyParticle(particle, tparticle);
		
		for(int i = 0; i < numVM; i++)
			for(int j = 0; j < numPM; j++)
			{
				if(!velocity[i][j])
				{
					if(Math.random() <= p1)
					{
						tparticle[i][j] = particle[i][j];
					}	
					else if(Math.random() <= p2 && Math.random() > p1)
					{
						tparticle[i][j] = lbest[i][j];
					}
					else if(Math.random() <= p3 && Math.random() > p2)
					{
						tparticle[i][j] = gbest[i][j];
					}
				}
			}
		if(this.isFeasible(this.matrixToArry(tparticle)))
		{
			double fitness = this.fitness(this.matrixToArry(tparticle));
			
			if(fitness > this.fitness[index])
			{
				this.fitness[index] = fitness;
				this.copyParticle(tparticle, lbest);
			}
			
			if(fitness > this.bestFitness)
			{
				this.bestFitness = fitness;
				this.copyParticle(tparticle, gbest);
			}
			this.copyParticle(tparticle, particle);
		}
	}
	
	private boolean[][] arrayToMatrix(int [] placement)
	{
		boolean [][] particle = new boolean [numVM][numPM];
		
		for(int i = 0; i < numVM; i++)
		{
			
			int pm = placement[i];
			if(pm != -1)
				particle[i][pm] = true;
		}
		
		return particle;
	}
	
	private int[] matrixToArry(boolean [][] particle)
	{
		int [] placement = new int[numVM];
		
		for(int i = 0 ; i < this.numVM; i++)
			for(int j = 0; j < this.numPM; j++)
			{
				if(particle[i][j])
					placement[i] = j;
			}
		
		return placement;
	}
	private void copyParticle(boolean[][] src,boolean [][]dst)
	{
		
			for(int i = 0 ; i < this.numVM; i++)
				for(int j = 0; j < this.numPM; j++)
					dst[i][j] = src[i][j];
	}// copyParticle
	
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
	
	
	
	private boolean generatePop(int[][] pop, int popSize) {
		
		int[][] useCPU = new int[popSize][numPM];
		int[][] useMEM = new int[popSize][numPM];
		boolean [] vmPlacement = new boolean[numVM];
		boolean first = true;
		boolean isLinkOverLoad = false;
		this.feasible = new boolean[popSize];
		int choose;
		int randomSize = popSize - popSize / 3;
		
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
				return true;
		
		return  false;
	}//generatePop
	
	private int nextPlacmented(boolean[] vmPlacement, int numVM2) {
		
		Random random = new Random();
		int choose = random.nextInt(numVM);
	
		while(vmPlacement[choose])
		{
			choose = ++choose % numVM;
		}
		
		return choose;
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
		
	}// end of updateLinkLoadByVM.
	
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
	}// end of isFeasible
	
	private void initParticlesPop(int[][]pop) {
		// TODO Auto-generated method stub
		this.fitness = new double[numParticle];
		this.bestIndividual = new int[numVM];
		Random  random = new Random();
	
		for(int k = 0; k < numParticle; k++)
			for(int i = 0; i < numVM; i++){
				for(int j = 0; j < numPM; j++){
					
					if(Math.random() > 0.5)
						v_velocity[k][i][j] = false;
					else
						v_velocity[k][i][j] = true;	
				
					}
				}
		
		double max_fitness = Double.MIN_VALUE;
		
	
		
		for(int k = 0; k < numParticle; k++)
		{	
			x_position[k] = this.arrayToMatrix(pop[k]);
			this.copyParticle(x_position[k], this.lBest[k]);
//			this.placeDetil("pso:", this.matrixToArry(this.lBest[k]));
			this.fitness[k] = this.fitness(pop[k]);
			
			if (this.fitness[k] > max_fitness)
			{	
				max_fitness = this.fitness[k];
				System.arraycopy(pop[k], 0, bestIndividual, 0, numVM);
				this.copyParticle(x_position[k], this.gBest);
				this.bestFitness = max_fitness;
			}	
		}
	}

	private double fitness(int[] c) {
		
		double fitness = 0 ;
		
		if(this.isFeasible(c))
		{
			fitness = this.SERVER_FULL_ENERGY / (this.totalEnergy(c) + sumTraffic(c));

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
			return fitness;
		}			
	}//end of fitness
	
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
	}//end of totalEnergy

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
	
	public double getFfTotalEnergy()
	{
		
		this.psoys_energy = this.totalEnergy(this.matrixToArry(gBest));
//		System.out.println(psoys_energy);
		return this.psoys_energy;
	}//end of getFfdTotalEnergy

}
