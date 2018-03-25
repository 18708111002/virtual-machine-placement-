package VMP_Algorithm;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;
public class IGAVP {
	

	private final static double SWITCH_IDLE_ENERGY = 147 ;

	
	private final static double SERVER2_CPU_MIPS = 24000;
	private final static double SERVER2_MEM = 16000;
	
	private final static double HIGH_SERVER_FULL_ENERGY = 550;//252;
	private final static double SERVER_FULL_ENERGY = 420;//252;
	private final static double SERVER_IDLE_ENERGY = SERVER_FULL_ENERGY * 0.7;
	private final static double HIGH_SERVER_IDLE_ENERGY = HIGH_SERVER_FULL_ENERGY * 0.7;
	
	
	private final static double PM_EDGE_CAPACITY = 1000 * 0.8;
	private final static double EDGE_AGGRE_CAPACITY = 1000 * 0.8;
	private final static double AGGRE_CORE_CAPACITY = 1000 * 0.8;
	
	private final static int TYPE_PORT_10MB = 10;
	private final static int TYPE_PORT_100MB = 100;
	private final static int TYPE_PORT_1000MB = 1000;
	
	private static double _10MB_ENERGY = 0.2;
	private static double _100MB_ENERGY = 0.4;
	private static double _1000MB_ENERGY = 1.1;
	
	private final double PR_MUTATION = 0.1;  
	private final double PR_CROSSOVER = 0.5;
	
	private final static double SERVER_CPU = 18080 * 2.5;
	private final static double SERVER_MEM = 64 * 1024;
	
	private int numVM;
	private int numPM;
	private int fatTreePod;
	
	private int core;
	private int aggre;
	private int edge;
	
	private boolean[] coreUsed; 
	
	private Pod[] pod;
	
	private double[] fitness;
	private int[] bestIndividual = null;

	private double worstFitness = Double.MAX_VALUE;
	private int worstIndex;
	private double bestFitness;
	private double totalEnergy;
	
	private int[] pmCPU; //the CPU capacity of pm.
	private int[] vmCPU; //the CPU requirement of vm.
	private int[] pmMEM; //the memory capacity of pm.
	private int[] vmMEM; //the memory requirement of vm.
	
	
	private double[][] c;

	private boolean[] feasible;
	private String fileName;
	private String file;
	private int[][] throughSwitch;
	
	private double maxFitness = Double.MIN_VALUE;

	public double getGaTotalEnergy()
	{
		return this.totalEnergy;
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
		
		this.coreUsed = new boolean[core];
		
		for(int i = 0; i < fatTreePod; i++)  // generate the fat-tree topology.
			this.pod[i] = new Pod(fatTreePod);
		
		this.core = (fatTreePod / 2) * (fatTreePod / 2);
		this.aggre = fatTreePod * fatTreePod / 2;
		this.aggre = fatTreePod * fatTreePod / 2;
	}
	
	public IGAVP(int numVm, int numPm) throws IOException
	{
		this.numPM = numPm;
		this.numVM = numVm;
		
		this.fileName = fileName;
		this.file = file;
		
		this.generateTopology(numPM);
		this.generateThroughSwitch();
		
		this.bestIndividual = new int[numVM];
			
	}
	
	public void getResource(int[] pmCPU, int[] pmMEM, int[] vmCPU, int[] vmMEM, double[][] c)
	{
		this.pmCPU = pmCPU;
		this.pmMEM = pmMEM;
		this.vmCPU = vmCPU;
		this.vmMEM = vmMEM;
		this.c = c;
	}
	
	private double weightArraySum(double [] weightArrays) {  
		double weightSum = 0;  
		for (double weightValue : weightArrays) {  
			weightSum += weightValue;  
		}  
		return weightSum;  
	 
	}  
	
	private void pairUp(int[] pairUp,int popSize)
	{
		WeightRandom random = new WeightRandom();
		double weightSum = weightArraySum(this.fitness);
		
		for(int i = 0; i < popSize; i++)   //pair up
		{
			pairUp[i] = random.getWeightRandom(this.fitness,weightSum);
			
			while(pairUp[i] == i)
				pairUp[i] = random.getWeightRandom(this.fitness,weightSum);
		}
	}
	
	private int getBestIndividual()
	{
		double f = Double.MIN_VALUE;
		int index = -1;
		
		for(int i = 0; i < this.fitness.length; i++)
		{
			if(this.fitness[i] > f)
			{
				f = this.fitness[i];
				index = i;
			}
		}
		
		return index;
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
	
	public int[] IGAVP(int popSize, int termination)
	{
		
		int[][]pop = new int[popSize][numVM];
		this.fitness = new double[popSize];
		
		this.generateIgavpPop(pop, popSize);
		
		
		FFA evmp = new FFA(numVM,numPM);
		evmp.getResource(getPmCPU(), getPmMEM(), getVmCPU(), getVmMEM(), getTraffic());
		pop[0] = evmp.FFA();
		
		this.fitness[0] = this.isgavpFitness(pop[0]);
		
		
		System.arraycopy(pop[0], 0,this.bestIndividual, 0, numVM);
		
		for(int i = 1 ; i < popSize; i++)
		{
			
			this.fitness[i] = this.isgavpFitness(pop[i]);	
	
			if(this.fitness[i] < this.worstFitness)
			{
				this.worstFitness = this.fitness[i];
				this.worstIndex = i;
			}
		}
		
		
		
		int count = 0;
		
		
		
		while(count < termination)
		{
			
			
			WeightRandom random = new WeightRandom();
			double weightSum = weightArraySum(this.fitness);
			
			int index1 = random.getWeightRandom(this.fitness, weightSum);
			int index2 = random.getWeightRandom(this.fitness, weightSum);
			
			int []newC = this.igavapCrossover(pop[index1],pop[index2]);
			this.igavapMutation(newC);
			
			for(int i = 0 ; i < popSize; i++)
			{
				this.fitness[i] = this.isgavpFitness(pop[i]);	
				
				if(this.fitness[i] < this.worstFitness)
				{
					this.worstFitness = this.fitness[i];
					this.worstIndex = i;
				}
			}
			
			System.arraycopy(newC, 0, pop[this.worstIndex], 0, numVM);
			count ++;
			
			
		}
		
		this.totalEnergy = this.totalEnergy(this.bestIndividual);
		
		return this.bestIndividual;
		
	}// END OF IGAVP
	
	private void igavapMutation(int [] c)
	{
		Random random =  new Random();
		
		for(int i = 0; i < c.length; i++)
		{
			if(Math.random() < this.PR_MUTATION)
				c[i] = random.nextInt(numPM);
		}		
	}
	
	private int [] igavapCrossover(int[]c1,int[]c2)
	{
		Random random = new Random();
		
		int[]tc1= new int[numVM];
		int[]tc2= new int[numVM];
		
		System.arraycopy(c1, 0, tc1, 0, numVM);
		System.arraycopy(c2, 0, tc2, 0, numVM);
		
		int postion = random.nextInt(numVM);
		
		for(int i = 0; i < postion; i++)
		{
			tc1[i] = tc2[i];
		}
		
		
		for(int i = postion; i < numVM; i++)
		{
			tc2[i] = tc1[i];
		}
		
		double f1 = this.isgavpFitness(tc1);
		double f2 = this.isgavpFitness(tc2);
		
		if(f1 > f2)
		{
			if(f1 > this.bestFitness)
				System.arraycopy(tc1, 0, this.bestIndividual, 0, numVM);
			
			return tc1;
		}
		else
		{
			if(f2 > this.bestFitness)
				System.arraycopy(tc1, 0, this.bestIndividual, 0, numVM);
			
			return tc2;
		}
		
	}// end of igavapCrossover
	
	private double  isgavpFitness(int[] c)
	{
		double fitness = 0;
		int count = 1;
		boolean [] usePM = new boolean[numPM];
		
		for(int i = 0; i < numPM;i++)
			usePM[i] = false;
		
		ArrayList<Integer> vmList = new ArrayList<Integer>();
		
		if(this.isFeasible(c))
		{
			
			
			for(int i = 0; i < this.numPM; i++)
			{
				int useCPU = 0;
				int useMEM = 0;
				
				for(int j = 0; j < this.numVM; j++)
				{
					usePM[c[j]] = true;
					
					if(c[j] == i)
					{
						vmList.add(j);
					}
				}
				for(int j = 0; j < vmList.size(); j++)
				{
					useCPU += this.vmCPU[j];
					useMEM += this.vmMEM[j];
				}
				
				fitness += (useCPU + useMEM); 	
			}
		}
		else
		{
			fitness = Double.MIN_VALUE;
		}
		
	
		for(int i = 0; i < numPM;i++)
		{
			if(usePM[i] == true)
			{
				count ++;
			}
		}
		
		return fitness / (count * 1.0);
	}
	
	public int[] GA(int popSize, int termination)
	{
		Random random = new Random();
		double curBestFitness = -1;
		int notImprove = 0;

		int[][]pop = new int[popSize][numVM];
		int[] pairUp = new int[popSize];   // pairUp[i] = j ;   i pair up j.
		this.fitness = new double[popSize];
			
		if(!this.generatePop(pop, popSize))
		{
//			System.out.println("GA no fesiable individual");
			return null;
		}
		
		this.initTopology();
		this.calculateFitness(popSize, pop);

		while(notImprove < termination)
		{
			this.initTopology();
			double averageFitness = 0;
			
			for(double fitness : this.fitness)
			{
				averageFitness += fitness;
			}
//			pair-up
			this.pairUp(pairUp, popSize);
			
			for(int i = 0; i < popSize; i++) //crossover
			{
				if(Math.random() < this.PR_CROSSOVER)
				{
					this.biasedUniformCrossover(this.fitness[i],fitness[pairUp[i]],pop[i], pop[pairUp[i]],
												pop[worstIndex],this.getWorstIndex());
				}
				
			}	
	
			for(int i = 0; i < popSize; i++)
			{
				this.mutation(pop[i],i);
			}
				
			int curBestIndex = this.getBestIndividual();
			curBestFitness = this.fitness[curBestIndex];
			
			
			
			if(curBestFitness > this.bestFitness)
			{
				System.arraycopy(pop[curBestIndex], 0, bestIndividual, 0, numVM);
				this.bestFitness = curBestFitness;
				notImprove = 0;
			}
			else if(curBestFitness != -1)
				notImprove ++;
			
//			System.out.println(this.bestFitness);
		}
		
		this.totalEnergy = this.totalEnergy(bestIndividual);
		return bestIndividual;
}

	private void mutation(int[] c,int index)
	{	
		Random random =  new Random();
		
		for(int i = 0; i < c.length; i++)
		{
			if(Math.random() < this.PR_MUTATION)
				c[i] = random.nextInt(numPM);
		}	
		
		this.correction(c);
		this.fitness[index] = this.fitness(c);
	}
	
	private void biasedUniformCrossover(double fit_i,double fit_j,int [] ci,int [] cj,int [] worst,int worstIndex)
	{
		double crossoverFitness = 0;
		int [] child = new int[numVM];
		for(int i = 0 ; i < ci.length; i++)
		{
			if(Math.random() < fit_i / (fit_i + fit_j))
			{
				child[i] = ci[i];
			}
			else
			{
				child[i] = cj[i];
			}
		}
		
//		this.correction(child);
		crossoverFitness = this.fitness(child);
		
		if(crossoverFitness > this.fitness[worstIndex])
		{
//			System.out.println("parent fitness:" +  fit_i + " " + fit_j);
//			System.out.println("best  fitness:" +  this.bestFitness );
//			System.out.println("crossover fitness:" + crossoverFitness);
//			System.out.println("worst fitness:" +  this.fitness[worstIndex] );
//			System.out.println();
			this.fitness[worstIndex] = crossoverFitness;
			System.arraycopy(child, 0, worst, 0, numVM);
		}

	}
	
	private void correction(int [] c)
	{
	
		ArrayList vm_unplacement =  new ArrayList();
		
		int [] useCPU = new int[numPM];
	    int [] useMEM = new int[numPM];
	  
	    for(int j = 0; j < numPM; j++)
	    {
	    	useCPU[j] = 0;
	    	useMEM[j] = 0;
	    }
		
		for(int i = 0; i < numVM ;i ++)
		{
			 if(c[i] == -1)
			 {
				 vm_unplacement.add(i);
				 continue;
			 }
			 
			 if(useCPU[c[i]] + this.vmCPU[i] > this.pmCPU[c[i]] ||
					 useMEM[c[i]] + this.vmMEM[i] > this.pmMEM[c[i]])
			 {		
				 vm_unplacement.add(i);
			 }
			 else 
			 {
				useCPU[c[i]] += this.vmCPU[i];
				useMEM[c[i]] += this.vmMEM[i];
			 }
		}
		
		for(Object vm : vm_unplacement)
		{
			int index = (int)vm;
			for(int i = 0; i < numPM ; i++)
			{
				if(useCPU[i] + this.vmCPU[index] <= this.pmCPU[i] &&
						 useMEM[i] + this.vmMEM[index] <= this.pmMEM[i])
				 {		
					 c[index] = i;
					 useCPU[i] += this.vmCPU[index];
					 useMEM[i] += this.vmMEM[index];
				 }
				
			}
		}
	
	}
	private void calculateFitness(int popSize, int[][] pop) {
		
		double max_fitness = Double.MIN_VALUE;
		double min_fitness = Double.MAX_VALUE;
		
		for(int i = 0; i < popSize; i++)
		{
			if(pop[i] == null)
				continue;
			
			this.fitness[i] = this.fitness(pop[i]);
			if (this.fitness[i] < min_fitness)
			{	
				min_fitness = this.fitness[i];
				this.worstIndex = i;
			}
			
			if (this.fitness[i] > max_fitness)
			{	
				max_fitness = this.fitness[i];
				System.arraycopy(pop[i], 0, bestIndividual, 0, numVM);
				this.bestFitness = max_fitness;
			}
		}
		
	}
	
	private int getWorstIndex()
	{
		int index = -1;
		double min_fitness = Double.MAX_VALUE;
		
		for(int i = 0; i < this.fitness.length; i++)
		{
			if (this.fitness[i] < min_fitness)
			{	
				min_fitness = this.fitness[i];
				index = i;
			}
		}
		
		return index;
		
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
	}
	
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
	
	private void generateThroughSwitch()
	{
		this.throughSwitch = new int[numPM][numPM];
		
		for(int i = 0; i < numPM; i++)
			for(int j = 0 ; j < numPM; j++)
				this.throughSwitch[i][j] = this.getThroughSwitch(i,j);
		
	}
	
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
	}
	
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
		
	}
	
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
	}
	

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
		
	}


	
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
	
	}
	
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
	
	private void updateLinkLoadByFFD(int index,int pm1,ArrayList vmRequire,int[] placement)
	{
		vmRequire vmr = (vmRequire) vmRequire.get(index);
		int vm = vmr.id;
		for(int i = 0; i < index; i++)
		{
			int vm2 = ((vmRequire)vmRequire.get(i)).id;
			double traffic = this.c[i][vm];
			
			int pm2 = placement[vm2];
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
			this.coreUsed[aggre1 * (this.fatTreePod / 2) + core1] = true;
			this.coreUsed[aggre2 * (this.fatTreePod / 2) + core2] = true;
			
			this.pod[pod1].pLinkE[pmIndex1][edge1] += traffic;
			this.pod[pod2].pLinkE[pmIndex2][edge2] += traffic;
			
			this.pod[pod1].eLinkA[edge1][aggre1] += traffic;
			this.pod[pod2].eLinkA[edge2][aggre2] += traffic;
			
			this.pod[pod1].eLinkA[aggre1][core1] += traffic;
			this.pod[pod2].eLinkA[aggre2][core2] += traffic;
					
		}
		
	}//
	
	private boolean generateIgavpPop(int[][] pop, int popSize)
	{
		int[][] useCPU = new int[popSize][numPM];
		int[][] useMEM = new int[popSize][numPM];
		boolean [] vmPlacement = new boolean[numVM];
		boolean first = true;
		boolean isLinkOverLoad = false;
		this.feasible = new boolean[popSize];
		int choose;
		int randomSize = popSize - 20;
		
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
		System.out.println("ok");
		
		
		for(int i = 0; i < popSize; i++)
			if(this.feasible[i])
				return true;
		
		return  false;
	}
	
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
		 
//		 for(int i = 0; i < numVM; i++)
//			 for(int j = 0; j < numPM; j++)
//			 {
//				 if(this.isLinkOverLoad(c[i], c[j], this.c[i][j]))
//					 return false;
//			 }
		 
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
	

	public void updateLinkLoadByVM(int vm, int pm1, int[] placement) {
		// TODO Auto-generated method stub
		
		for(int i = 0; i < vm; i++)
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
			
			this.coreUsed[aggre1 * (this.fatTreePod / 2) + core1] = true;
			this.coreUsed[aggre2 * (this.fatTreePod / 2) + core2] = true;
			
			
			this.pod[pod1].aLinkC[aggre1][core1] += traffic;
			this.pod[pod2].aLinkC[aggre2][core2] += traffic;
					
		}
		
	}// end of updateLinkLoadByVM.
	
	private int getPmIndex(int podNum,int edge,int pm)
	{
		return ( pm - (podNum * fatTreePod * fatTreePod / 4 + edge * fatTreePod / 2));
		
	}//end of getPmIndex
	
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
		
		boolean p = false;
		
//		System.out.println(podNum);
		Pod pod = this.pod[podNum];
		
		for(int i = 0; i < pod.eLinkA[edge].length;i++)
		{
			if(pod.eLinkA[edge][i] < minTraffic)
			{
				minTraffic = pod.eLinkA[edge][i];
				aggre = i;
			}
		}
		
		return aggre;
	}//end of getAggre
	
	private int getCore(int podNum,int aggre)
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
		
		this.coreUsed[aggre * (this.fatTreePod / 2) + core] = true;
		
		return core;
	}//end of getCore
	
	private int nextPlacmented(boolean[] vmPlacement, int numVM2) {
		
		Random random = new Random();
		int choose = random.nextInt(numVM);
	
		while(vmPlacement[choose])
		{
			choose = ++choose % numVM;
		}
		
		return choose;
	}
}
