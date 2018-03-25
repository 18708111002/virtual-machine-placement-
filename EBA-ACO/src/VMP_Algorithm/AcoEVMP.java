package VMP_Algorithm;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

public class AcoEVMP {
	
	public double aco_bw = 0;
	public int[] ffdPlacement;
	
	private int[] bestPlacement;
	private List<Integer> vmList = new ArrayList<Integer>();
 
	private int numVM;
	private int numPM;
	
	
	private double [][] c;
	private int[][] throughSwitch;
	public double vmTotalTraffic = 0;
	private int [] vmCPU;
	private int [] vmMEM;
	
	private int [] pmCPU;
	private int [] pmMEM;
	
	private double totalEnergy = Double.MAX_VALUE;
	private double curBestEnergy;

	private Object fileName;
	private String file;

	private int fatTreePod;
	public int core;
	public int aggre;
	public int edge;
	private Pod[] pod;

	private boolean[] coreUsed;
	private double maxTraffic = 0;
	private double GAMMA = 1;

	
	private final static double SERVER1_CPU_MIPS = 40000;
	private final static double SERVER1_MEM = 32000;
	
	private final static double SERVER2_CPU_MIPS = 24000;
	private final static double SERVER2_MEM = 16000;
	
	private final static double HIGH_SERVER_FULL_ENERGY = 550;//252;
	private final static double SERVER_FULL_ENERGY = 420;//252;
	private final static double SERVER_IDLE_ENERGY = SERVER_FULL_ENERGY * 0.7;
	private final static double HIGH_SERVER_IDLE_ENERGY = HIGH_SERVER_FULL_ENERGY * 0.7;
	
	
	private final static double SWITCH_IDLE_ENERGY = 147 ;
	
//	private final static double CORE_SWITCH_FULL_ENERGY = 650 ;
	private final static double CORE_SWITCH_IDLE_ENERGY = 555 ;
	private final static double EDGE_SWITCH_IDLE_ENERGY = 150 ;
	
	
	private final static int TYPE_PORT_10MB = 10;
	private final static int TYPE_PORT_100MB = 100;
	private final static int TYPE_PORT_1000MB = 1000;
	
	private static double _10MB_ENERGY = 0.2;
	private static double _100MB_ENERGY = 0.4;
	private static double _1000MB_ENERGY = 1.1;
	
	private final static double PM_EDGE_CAPACITY = 10000 * 0.8;
	private final static double EDGE_AGGRE_CAPACITY = 1000 * 0.8;
	private final static double AGGRE_CORE_CAPACITY = 1000 * 0.8;
	private final static double COMMUNICATION_BUFFER =  50;


	
	private final static double INIT_PHEROMONE = 0.0001;
	private final static double INIT_HEURISTIC = 0;
	
	
	public  final static double BETA = 2 ;
	public  final static double ALPHA = 1;
	private final static double EPSILON = 0.0001;
	private final static double RHO_LOCAL = 0.4;
	private final static double RHO_GLOBAL = 0.35;
	private final static double VM_PROPORTION = 0.3; 
	
	private final static int    BIG_FLOW_NUM = 10;
	
	private final double PR_BIAS = 0.1; 
	
	public double switchChassis = 0;
	private double maxAntEnergy;


	public AcoEVMP(int numVm, int numPm,String fileName,String file) throws IOException
	{
		this.numPM = numPm;
		this.numVM = numVm;
		this.bestPlacement = new int[numVM];

		this.fileName = fileName;
		this.file = file;
		
		this.getVmRequire();
		this.getPmCapacity();
		this.getTrafficRequire();
//		
		this.generateTopology(numPm);	
		this.generateVmList();
		this.generateThroughSwitch();
		
	}// end of ACO_EVMP
	

	private void generateThroughSwitch()
	{
		this.throughSwitch = new int[numPM][numPM];
		
		for(int i = 0; i < numPM; i++)
			for(int j = 0 ; j < numPM; j++)
				this.throughSwitch[i][j] = this.getThroughSwitch(i,j);
		
	}
	private void generateVmList()
	{
		List<Integer> vmList = new ArrayList<Integer>();
		
		ArrayList vmTraffic = new ArrayList();
		boolean[] vmAdded = new boolean [numVM];
		
		for(int i = 0; i < numVM; i++)
			vmAdded[i] = false;
		
		for(int i = 0; i < numVM; i++)
			for(int j = i + 1; j < numVM; j++)
			{
				vmTraffic.add(new VMtraffic(i,j,this.c[i][j]));
			}

		vmTraffic.sort(new SortByTraffic());
		
		for(Object vt:vmTraffic)
		{
			int vm1 = ((VMtraffic)vt).vm1;
			int vm2 = ((VMtraffic)vt).vm2;
			
			if(!vmAdded[vm1]) 
			{
				this.vmList.add(vm1);
				vmAdded[vm1] = true;
			}
			
			if(!vmAdded[vm2]) 
			{
				this.vmList.add(vm2);
				vmAdded[vm2] = true;
			}	
		}
		
		for(int i = 0; i < numVM;i++)
		{
			if(!vmAdded[i])
			{
				this.vmList.add(i);
				vmAdded[i] = true;
			}
		}	
	}
	
	private void getTrafficRequire() throws IOException
	{
		this.c = new double[numVM][numVM];	
		RequirementWR.readRequirement(c, fileName + "TRAFFIC\\" + file);
		
		 for (int x = 0; x < this.c.length; x++) {  
	            for (int y = 0; y < this.c[x].length; y++) {  
	            	
	            	if(this.c[x][y] > this.maxTraffic)
	            		this.maxTraffic = this.c[x][y];
	            	
	                this.vmTotalTraffic += this.c[x][y] ;  
	            }  
		 }
	}
	
	private void getVmRequire() throws IOException
	{
		this.vmCPU = new int[numVM];
		this.vmMEM = new int[numVM];
		
		RequirementWR.readRequirement(vmCPU, fileName + "VM-CPU\\" + file);
		RequirementWR.readRequirement(vmMEM, fileName + "VM-MEM\\" + file);
	}
	
	private void getPmCapacity() throws IOException
	{
		this.pmCPU = new int[numPM];
		this.pmMEM = new int[numPM];
		
		
		RequirementWR.readRequirement(pmCPU,fileName + "PM-CPU\\" + file);
		RequirementWR.readRequirement(pmMEM,fileName + "PM-MEM\\" + file);
				
	}
	private void generateVmRequire()
	{
		
		Random random = new Random();
		this.vmCPU = new int[numVM];
		this.vmMEM = new int[numVM];
		
		
		for(int i = 0; i < numVM; i++)
		{
			vmCPU[i] = 1500 + random.nextInt(3500 + 1);
//			vmCPU[i] = (int) (3500 + Math.random() * 1000);
			vmMEM[i] = 1500 + random.nextInt(3500 + 1);
//			vmMEM[i] = (int) (4500 + Math.random() * 1000);
		}
	}// end of generateVmRequire
	private void generatePmCapacity()
	{
		this.pmCPU = new int[numPM];
		this.pmMEM = new int[numPM];
		
		
		boolean[] pmAlloc = new boolean[this.numPM];
		
		for(int i = 0 ; i < this.numPM ; i++)
			pmAlloc[i] = false;
		
		int pmCount = (int) (this.numPM * 0.3);
		int anotherPmCount = this.numPM - pmCount;
		
		Random random = new Random();
		
		while(pmCount > 0)
		{
			int indexPm = random.nextInt(this.numPM);
			
			while(pmAlloc[indexPm])
				indexPm = random.nextInt(this.numPM);
			
			pmCPU[indexPm] = (int) this.SERVER1_CPU_MIPS;
			pmMEM[indexPm] = (int) this.SERVER1_MEM;
			
			pmAlloc[indexPm] = true;
			pmCount --;
				
		}
		
		
		while(anotherPmCount > 0)
		{
			int indexPm = random.nextInt(this.numPM);
			
			while(pmAlloc[indexPm])
				indexPm = random.nextInt(this.numPM);
			
			pmCPU[indexPm] = (int) (this.SERVER2_CPU_MIPS);
			pmMEM[indexPm] = (int) (this.SERVER2_MEM);
			
			pmAlloc[indexPm] = true;
			anotherPmCount --;
				
		}
	}// end of generatePmCapacity
	private void generateVmTraffic()
	{
		java.text.DecimalFormat   df   =new   java.text.DecimalFormat("#.00");  
		Random random = new Random();
		List<Integer> vmList = new ArrayList<Integer>();
		
		int count = numVM;
		double traffic;

		this.c = new double[numVM][numVM];
		
		for(int i = 0; i < numVM; i++)
			vmList.add(i);
		
		while(vmList.size() > 0)
		{
			int num = 3 + random.nextInt(10);
			int index = random.nextInt(vmList.size());
			int vm1 = vmList.get(index);
			vmList.remove(index);
		
			
			for(int j = 0; j < num; j++)
			{
				if(vmList.size() > 0)
				{
					int indexK = random.nextInt(vmList.size());
					int k = vmList.get(indexK);
					vmList.remove(indexK);
					
					if(vmList.size() > BIG_FLOW_NUM )
					{
						traffic = 1 +  Math.random() * 5;
			
						this.c[vm1][k] = Double.valueOf(df.format(traffic));
						this.c[k][vm1] = this.c[vm1][k];
					}
					else
					{
						traffic = 10 + Math.random() * 100;
		
						this.c[vm1][k] = Double.valueOf(df.format(traffic));
						this.c[k][vm1] = this.c[vm1][k];
						
						if(this.maxTraffic < this.c[k][vm1])
							this.maxTraffic = this.c[k][vm1];
					}
				}
			}
			
		}
		
		for (int x = 0; x < this.c.length; x++) {  
            for (int y = 0; y < this.c[x].length; y++) {  
            	
                this.vmTotalTraffic += this.c[x][y] ;  
            }  
	 }
		
	}// end of generateVmTraffic
		

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
		
		

		
			
		this.switchChassis = this.core * this.CORE_SWITCH_IDLE_ENERGY + 
							 (this.aggre + this.edge) * this.EDGE_SWITCH_IDLE_ENERGY;
		
	}
	public int [] ACO(HashMap<String , Object> parameters )
	{	
		int ants = (int)parameters.get("ants");
		int maxGen = (int)parameters.get("maxGen");
		double q = (double)parameters.get("q");
		
		
		double[] perAntEnergy = new double[ants];

		int [][] placeSet = new int [ants][numVM];
		
		int [][] localBestPlaceSet = new int[ants][numVM];
		double[] perLocalAntEnergy = new double[ants];
		
		for(int i = 0; i < ants; i++)
			perLocalAntEnergy[i] = Double.MAX_VALUE;
		
		
		double [][] pheromone = new double [numVM][numPM];
		double [][] heuristic = new double [numVM][numPM];
		
		int[] useCPU = new int[numPM];
		int[] useMEM = new int[numPM];
		boolean[] usePM  = new boolean[numPM];
		int[] bestPlacement;
		
//		int[] guidePlacement = this.generateGuidePlacement(ants);
		this.initTopology();
		this.initACO(pheromone,heuristic); // init parameter.
		int curGen = 0;
		int pm;
		
		while(curGen < maxGen)
		{
			for(double energy : perAntEnergy) 
				energy = 0;
			
			this.maxAntEnergy = 0 ;
			
			for(int k = 0; k < ants; k++ )
			{
				// init the server's resoures used by previous ant.
				this.initServerUsedResource(useCPU,useMEM,usePM,placeSet[k]);
				
				for(int i = 0; i < this.numVM; i++)
				{
					ArrayList serverList = this.getServerList(i,useCPU,useMEM,placeSet[k]); // servers can hold the current VM i.
					
					if(serverList.size() > 0)
					{
						this.updateHeuristic(serverList,i,heuristic,useCPU,useMEM,usePM,placeSet[k]);
						
						if(Math.random() < q)
							pm = this.selectServerByExp(placeSet[k],i,serverList,pheromone,heuristic,usePM,useCPU,useMEM);
						else
							pm = this.selectServerByPr(placeSet[k],i,serverList,pheromone,heuristic,usePM,useCPU,useMEM);
					
						this.updateLocalPheromone(pheromone,i,pm);
						this.updateLinkLoad(i,pm,placeSet[k]);
							
					}
					else
					{
						System.out.println("no hold server");
						return null;
					}
	
				} // end of VM
				this.calculateTotalEnergy(placeSet[k],perAntEnergy,k);
				
//				if(perAntEnergy[k] < perLocalAntEnergy[k])
//				{
//					perLocalAntEnergy[k] = perAntEnergy[k];
//					System.arraycopy(placeSet[k], 0, localBestPlaceSet[k], 0, numVM);
//				}
				
				
			} // end of ants
			
			bestPlacement = this.getBestPlacement(placeSet,perAntEnergy);
			this.updateGlobalPheromone(bestPlacement,pheromone,curBestEnergy);
			
//			double best
			double[] antsFitness = new double[ants];
			int bestIndex = this.calculateAntsFitness(ants, placeSet, antsFitness);
			
			int j = bestIndex;
			for(int i = 0; i < ants ; i++)
//				for(int j = 0; j < ants; j++)
				{
					if(i != j)
					{
						this.antsCommunication(antsFitness[i], antsFitness[j], placeSet[i], placeSet[j], i,perAntEnergy);
					}
				}
			
			bestPlacement = this.getBestPlacement(placeSet,perAntEnergy);
			this.updateGlobalPheromone(bestPlacement,pheromone,curBestEnergy);
			
//		    bestPlacement = this.getBestPlacement(localBestPlaceSet,perLocalAntEnergy);
//			this.updateGlobalPheromone(bestPlacement,pheromone,curBestEnergy);
//			
//			double[] antsFitness = new double[ants];
//			this.calculateAntsFitness(ants, localBestPlaceSet, antsFitness);
//			
//			for(int i = 0; i < ants ; i++)
//				for(int j = 0; j < ants; j++)
//				{
//					if(i != j)
//					{
//						this.antsCommunication(antsFitness[i], antsFitness[j], localBestPlaceSet[i], localBestPlaceSet[j], i,perLocalAntEnergy);
//					}
//				}
//			
//			bestPlacement = this.getBestPlacement(localBestPlaceSet,perLocalAntEnergy);
//			this.updateGlobalPheromone(bestPlacement,pheromone,curBestEnergy);
			
			curGen ++;
		} // end of while
		
		return this.bestPlacement;
	}// end of ACO
	
	
	private void antsCommunication(Double fit_i,Double fit_j,int [] ant_i,int [] ant_j,int ant, double[] perAntEnergy)
	{
		
		double newSolutionFitness = 0;
		int [] newSolution  = new int[numVM];
		for(int i = 0 ; i < ant_i.length; i++)
		{
			if(Math.random() < fit_i / (fit_i + fit_j))
			{
				newSolution[i] = ant_i[i];
			}
			else
			{
				newSolution[i] = ant_j[i];
			}
		}
		
		this.correction(newSolution);
		
		this.initTopology();
		for(int i = 0; i < numVM; i++)
		{
			int vm = this.vmList.get(i);
			int pm1 = newSolution[vm];
			this.updateLinkLoad(i, pm1, newSolution);
		}// update link load.
		
		if(this.isFeasible(newSolution))
		{
			newSolutionFitness = this.fitness(newSolution);

			if(newSolutionFitness > fit_i)
			{
				System.arraycopy(newSolution, 0, ant_i, 0, numVM);
				perAntEnergy[ant] = this.totalPmEnergy(newSolution);
				fit_i = newSolutionFitness;
			}
		}
		else
		{
//			System.out.println("infeasible");
		}
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
	}// end of correction.
	
	
	
	private int calculateAntsFitness(int ants, int[][] antsSolution,double [] antsFitness) {
		
	
		
		
		double max_fitness = Double.MIN_VALUE;
		int maxIndex = 0;
		
		
		for(int i = 0; i < ants; i++)
		{
			
//			System.out.println(this.getMeanMaxServerEnergy(antsSolution[i]));
//			antsFitness[i] = ( 1.0 ) / (this.totalPmEnergy(antsSolution[i]) + sumTraffic(antsSolution[i]));
			antsFitness[i] = ( 1.0 ) / (this.totalPmEnergy(antsSolution[i]));
			
			
			if (antsFitness[i] > max_fitness)
			{	
				max_fitness = antsFitness[i];
				maxIndex = i;
				
			}
		}
		
		
		return maxIndex;
		
		
		
		
//		double minVlaue = Double.MAX_VALUE;
//		int minIndex = 0;
////		
////		
//		for(int i = 0; i < ants; i++)
//		{
//		
//			antsFitness[i] = this.totalPmEnergy(antsSolution[i]) + sumTraffic(antsSolution[i]);
//			
//			
//			if (antsFitness[i] < minVlaue)
//			{	
//				minVlaue = antsFitness[i];
//				minIndex = i;
//				
//			}
//		}
////		
////		
//		return minIndex;
		
		
		
	}// end of calculateAntsFitness.
	
	private int activePm(int [] placement)
	{
		int total = 0;
		
		boolean [] pmEmpty = new boolean[numPM];
		
		
		for(int i = 0; i < numPM; i++)
		{
			pmEmpty[i] = true;
		}
		
		for(int i = 0; i < numVM; i++)
		{
			pmEmpty[placement[i]] = false;
		}
		
			
		for(int i = 0; i < numPM; i++)
		{
			if(!pmEmpty[i])
			{
				total ++;
			}
		}
		
		
		return total;
		
	}
	
	private int[] generateGuidePlacement(int ants) {
		// TODO Auto-generated method stub
		
		int [][]guidePlacement = new int[ants][numVM];
		boolean [] correctPlacemen = new boolean[ants];
		int[][] useCPU = new int[ants][numPM];
		int[][] useMEM = new int[ants][numPM];
		boolean [] vmPlacement = new boolean[numVM];
		boolean first = true;
		boolean isLinkOverLoad = false;
		int choose;
		
		for(int k = 0;  k < numVM; k++ )
		{
			vmPlacement[k] = false;	
		}
		
		
		for(int i = 0; i < ants; i++)
			for(int j = 0; j < numVM; j++)
				guidePlacement[i][j] = -1;
		
		for(int i = 0; i < ants; i++)
			for(int j = 0; j < numPM; j++)
			{
				useCPU[i][j] = 0;
				useMEM[i][j] = 0;
			}
	
		for(int i = 0; i < ants; i++)
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
							int pm1 = guidePlacement[i][t];
							
							if(this.isLinkOverLoad(k, pm1, this.c[j][vm]))
							{
//								System.out.println("guideplacement link is not enough");
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
						guidePlacement[i][choose] = k;
						vmPlacement[choose] = true;
						first = false;		
						
						this.updateLinkLoadByGuide(choose, k, guidePlacement[i], vmPlacement);
					}
				}					
			}	
			for(int t = 0; t < numVM; t++)
			{				
				vmPlacement[t] = false;	
			}
			first = true;				
		}
		
		double maxFitness = Double.MIN_VALUE;
		int maxIndex = -1;
		boolean allPlacement;
		
		for(int i = 0; i < ants; i++)
		{
			allPlacement = true;
			
			for(int j = 0; j < numVM; j++)
			{
				if(guidePlacement[i][j] == -1)
				{
					allPlacement = false;		
					break;
				}
			}
			if(allPlacement)
			{
				double fitness = this.fitness(guidePlacement[i]);
				if(maxFitness < fitness)
				{
					maxFitness = fitness;
					maxIndex = i;
				}
			}
				
		}
		
		if(maxIndex == -1)
			return null;
		
		return guidePlacement[maxIndex];
	}// end of generateGuidePlacement
	
	private void updateLinkLoadByGuide(int vm, int pm1, int[] placement, boolean[] vmPlacement) {
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
		}
		
	}

	private double fitness(int[] placement)
	{
		double pmEnergy = this.totalPmEnergy(placement);
//		double switchEnergy = this.totalSwitchEnergy(placement);
		double sumTraffic = sumTraffic(placement);
		double fitness = this.SERVER_FULL_ENERGY / (pmEnergy + sumTraffic);
		
//		System.out.println("pm:" + pmEnergy);
//		System.out.println("switch:" + switchEnergy);
//		System.out.println("traffic:" + sumTraffic);
		
		return fitness;
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
			
	}// end of sumTraffic
	
	private int nextPlacmented(boolean[] vmPlacement, int numVM2) {
		
		Random random = new Random();
		int choose = random.nextInt(numVM);
	
		while(vmPlacement[choose])
		{
			choose = ++choose % numVM;
		}
		
		return choose;
	}//end of nextPlacmented
	private void calculateTotalEnergy(int[] placement, double[] perAntEnergy, int k) 
	{
		double pmEnergy = this.totalPmEnergy(placement);
		double switchEnergy = this.totalSwitchEnergy(placement);
		
	
//		System.out.println("ACO"  +  pmEnergy);
//		System.out.println("ACO"  +  switchEnergy);
		
		perAntEnergy[k] =  pmEnergy + switchEnergy;
		
		if(perAntEnergy[k] == 0)
			perAntEnergy[k] = Double.MAX_VALUE;
		
		if(this.maxAntEnergy < perAntEnergy[k])
			this.maxAntEnergy = perAntEnergy[k];
		
	}//end of calculateTotalEnergy

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

	
	
	private void updateLinkLoad(int index, int pm1, int[] placement) {

		int vm = this.vmList.get(index);
		for(int i = 0; i < index; i++)
		{
			int vm2 = this.vmList.get(i);
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
			
			
	
			this.coreUsed[aggre1 * (this.fatTreePod / 2) + core1] = true;
			this.coreUsed[aggre2 * (this.fatTreePod / 2) + core2] = true;
			
			
			this.pod[pod1].aLinkC[aggre1][core1] += traffic;
			this.pod[pod2].aLinkC[aggre2][core2] += traffic;
					
		}
		
	}

	private void initACO(double[][] pheromone, double[][] heuristic)
	{
		for(int i = 0; i < numVM; i++)
			for(int j = 0; j < numPM; j++)
			{
				pheromone[i][j] = this.INIT_PHEROMONE;
				heuristic[i][j] = this.INIT_HEURISTIC;
			}
		
//		if(guidePlacement != null)
//			this.initPheromone(guidePlacement, pheromone);
		
	}// end of initACO	
	
	private void initPheromone(int[] placement, double[][] pheromone)
	{
		
		this.curBestEnergy = this.totalPmEnergy(placement) + this.totalSwitchEnergy(placement);
		ArrayList [] serverPlacement =  new ArrayList[numPM];
		int updateNum = 0;
		RankIndex rankIndex = new RankIndex();
		ArrayList<Double> lst = new ArrayList<Double>();
		double[] fitnessArry = new double[numPM];
		boolean [] updatedPM = new boolean[numPM];
		
		for(int k = 0;  k < numPM; k++)
		{
			serverPlacement[k] = new ArrayList();
			updatedPM[k] = false;	
		}
		
		for(int j = 0; j < numVM ;j ++)
		{
			serverPlacement[placement[j]].add(j);	
		}
		
		for(int i = 0; i < numPM; i++)
		{
			double fitness = this.evaluateSinglePM(serverPlacement[i], i);
		
			lst.add(fitness);
			fitnessArry[i] = fitness;
		
			if(fitness < (0.03 + 0.035))
				updateNum ++;	
//			System.out.println(this.evaluateSinglePM(serverPlacement[i], i) + "pm " + i + "updateNum" + updateNum);
		}
	
		if(updateNum > 0)
		{
			int [] index = new int[updateNum];
			
			for(int i = 0; i < index.length; i++)
	        	index[i] = i;
			
			rankIndex.minN(lst, updateNum, index);
			
			
			for(int i : index)
			{
				updatedPM[i] = true;
				for(Object vm: serverPlacement[i])
				{
					
					pheromone[(int)vm][i] = (1 - this.RHO_GLOBAL) * pheromone[(int)vm][i] + this.RHO_GLOBAL* (10000 / curBestEnergy + this.INIT_PHEROMONE / fitnessArry[i]);
//					System.out.println(pheromone[(int)vm][i]);
				}		
			}	
		}//end of if		
	}// end of initPheromone
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
	private void initServerUsedResource(int[] useCPU, int[] useMEM, boolean[] usePM, int[] placement)
	{
		
		this.initTopology();
		
		for(int i = 0; i < numPM; i++)
		{
			useCPU[i] = 0;
			useMEM[i] = 0;
			usePM[i]  = false;
		}
		
		for(int j = 0; j < numVM; j++)
			placement[j] = -1;
		
	}// end of initServerUsedResource
	private void updatePmResource(int vm, int pm, int[] useCPU, int[] useMEM)
	{
		useCPU[pm] += this.vmCPU[vm];
		useMEM[pm] += this.vmMEM[vm];
		
	}// end of updatePmResource
	private ArrayList getServerList(int index, int[] useCPU, int[] useMEM, int[] placement) // return servers can hold the VM.
	{
		int vm = this.vmList.get(index);
		ArrayList serverList = new ArrayList();
		
		for(int j = 0; j < numPM; j++)
		{
			if(	this.isResouceSatisfy(vm,j,useCPU,useMEM) &&
				this.isNetworkSatisfy(index,j,placement)	)
			{
				serverList.add(j);
			}
				
		}
		
		return serverList;	
	}// end of getServerList

	private int selectServerByExp(int[] placement, int index, ArrayList serverList, double[][] pheromone, double[][] heuristic, boolean[] usePM, int[] useCPU, int[] useMEM)
	{
		double maxExp = Double.MIN_VALUE;
		int pm = -1;
		int vm = this.vmList.get(index);
		int maxIndex = 0;
		double []exp = new double[serverList.size()];
		List<Integer> candidateServerIndex = new ArrayList<Integer>();
	
		
		for(int i = 0; i < serverList.size(); i++ )
		{	
			pm = (int) serverList.get(i);		
//			if(pheromone[vm][pm] != this.INIT_PHEROMONE)
//			{
//				System.out.println("pheromone:" +pheromone[vm][pm] * 10 );
//				System.out.println("heuristicOne:" + heuristic[vm][pm] * 10);
//			}
			exp[i] = Math.pow(pheromone[vm][pm] * 100 ,ALPHA) * Math.pow(heuristic[vm][pm] * 100,BETA);
//			System.out.println("pheromone:" + Math.pow(pheromone[vm][pm]  ,ALPHA));
//			System.out.println("heuristic:" + Math.pow(heuristic[vm][pm] * 100,BETA));
//			System.out.println(exp[i]);
			if(maxExp < exp[i])
			{
				placement[vm] = (int) serverList.get(i);	
				maxExp = exp[i];
				maxIndex = i;
			}
		}
		
		for(int i = 0; i < serverList.size(); i++)
		{	
			pm = (int) serverList.get(i);
			
			if(exp[i] > 0.7 * maxExp)
			{
				candidateServerIndex.add(i);
			}
		}
			
//		maxExp =  maxExp + this.heuristicFour(placement);
		maxExp =   maxExp * this.heuristicThree((int) serverList.get(maxIndex), vm, maxExp);
		
		for(int i : candidateServerIndex)
		{
//			placement[vm] = (int) serverList.get(i);
//			exp[i] = exp[i] + this.heuristicFour(placement);
//			placement[vm] = (int) serverList.get(maxIndex);
			
			exp[i] = exp[i] * this.heuristicThree((int) serverList.get(i), vm, exp[i]);
			
			if(maxExp < exp[i])
			{
//				System.out.println("placement");
				placement[vm] = (int) serverList.get(i);
				maxExp = exp[i];
				maxIndex = i;
			}
		}
		
		this.updatePmResource(vm,placement[vm],useCPU,useMEM);
		
		return placement[vm];
		
	}// end of selectServerByExp
	private int selectServerByPr(int[] placement, int index, ArrayList serverList, double[][] pheromone, double[][] heuristic, boolean[] usePM, int[] useCPU, int[] useMEM)
	{
		double []pr = new double [serverList.size()];
		int vm = this.vmList.get(index);
		
		for(int i = 0; i < serverList.size(); i++)
		{
			int pm = (int)serverList.get(i);
			
//			if(pheromone[vm][pm] != this.INIT_PHEROMONE)
//			{
//				System.out.println("pheromone:" +pheromone[vm][pm] * 10 );
//				System.out.println("heuristicOne:" + heuristic[vm][pm] * 10);
//				System.out.println(Math.pow(pheromone[vm][pm] * 10 ,ALPHA) *  Math.pow(heuristic[vm][pm] * 10 ,BETA));
//				System.out.println();
//				System.out.println();
//				
//			}
//			
			pr[i] = Math.pow(pheromone[vm][pm] * 100 ,ALPHA) *  Math.pow(heuristic[vm][pm] * 100 ,BETA);
//			pr[i] = Math.pow(pheromone[vm][pm] ,ALPHA);			
		}
		
		WeightRandom random = new WeightRandom();
		double weightSum = weightArraySum(pr);
		
		int pm = (int) serverList.get(random.getWeightRandom(pr, weightSum));
		placement[vm] = pm;
		usePM[pm] = true;
		
//		System.out.println(placement[vm]);
		
		this.updatePmResource(vm,placement[vm],useCPU,useMEM);
		
		return pm;
		
	}// end of selectServerByPr
	private void updateLocalPheromone(double[][] pheromone, int index, int pm)
	{
		int vm = this.vmList.get(index);
		pheromone[vm][pm] = (1 - this.RHO_LOCAL) * pheromone[vm][pm] + this.RHO_LOCAL * this.INIT_PHEROMONE;	
		
	}// end of updateLocalPheromone
	
	private int[] getBestPlacement(int[][] placeSet, double[] perAntEnergy)
	{
		int[] bestPlacement = new int[numVM];
		double minEnergy = Double.MAX_VALUE;
		double minUseSwitch = Double.MAX_VALUE;
		double curEnergy;
		double curUseSwitch;
		double minFitness = Double.MAX_VALUE;
		int switches = core + edge + aggre;
		
		for(int k = 0; k < placeSet.length; k++)
		{  
			int [] placement = placeSet[k];
			curEnergy = perAntEnergy[k];
		
			
//			double fitness = curEnergy *(1 + ((curUseSwitch  - minUseSwitch) / minUseSwitch - 0.1) / 10);
//			
//			if(minFitness > fitness)
//			{
//				minFitness = fitness;
//				minEnergy = curEnergy;
//				minUseSwitch = curUseSwitch;
//				System.arraycopy(placement, 0, bestPlacement, 0, numVM);
//			}
//	
			if(minEnergy > curEnergy )
			{
				minEnergy = curEnergy;
				System.arraycopy(placement, 0, bestPlacement, 0, numVM);
			}
			
			
//			if(minEnergy > curEnergy && (curUseSwitch  - minUseSwitch) / minUseSwitch < 0.1)
//			{
//				minEnergy = curEnergy;
//				minUseSwitch = curUseSwitch;
//				System.arraycopy(placement, 0, bestPlacement, 0, numVM);
//			}
		}	
		
		if(minEnergy < this.totalEnergy)
		{
			this.totalEnergy = minEnergy;
			System.arraycopy(bestPlacement, 0, this.bestPlacement, 0, numVM);
		}
		
		curBestEnergy = minEnergy;
		
		return bestPlacement;	
	}// end of getBestPlacement
	private void updateGlobalPheromone(int[] placement, double[][] pheromone, double curBestEnergy)
	{
		
//		ArrayList [] serverPlacement =  new ArrayList[numPM];
//		int updateNum = 0;
//		RankIndex rankIndex = new RankIndex();
//		ArrayList<Double> lst = new ArrayList<Double>();
//		double[] fitnessArry = new double[numPM];
//		boolean [] updatedPM = new boolean[numPM];
		
//		for(int k = 0;  k < numPM; k++)
//		{
//			serverPlacement[k] = new ArrayList();
//			updatedPM[k] = false;	
//		}
		
//		for(int j = 0; j < numVM ;j ++)
//		{
//			serverPlacement[placement[j]].add(j);	
//		}
//		
//		int count_2 = 0;
//		int count_1 = 0;
//		
//		for(int i = 0; i < numPM; i++)
//		{
//			double fitness = this.evaluateSinglePM(serverPlacement[i], i);
//		
//			lst.add(fitness);
//			fitnessArry[i] = fitness;
//		
//			if(fitness < (0.03 + 0.035))
//				updateNum ++;	
////			System.out.println(this.evaluateSinglePM(serverPlacement[i], i) + "pm " + i + "updateNum" + updateNum);
//		}
//	
//		if(updateNum > 0)
//		{
//			int [] index = new int[updateNum];
//			
//			for(int i = 0; i < index.length; i++)
//	        	index[i] = i;
//			
//			rankIndex.minN(lst, updateNum, index);
//			
//			
//			for(int i : index)
//			{ 
//				for(Object vm: serverPlacement[i])
//				{
//				
//					if(fitnessArry[i] > 0) 
//					{			
//	//					double bfPhr = pheromone[(int)vm][i];
//	//					System.out.println(this.INIT_PHEROMONE / fitnessArry[i] + " fitness"+ fitnessArry[i]);
//						pheromone[(int)vm][i] = (1 - this.RHO_GLOBAL) * pheromone[(int)vm][i] + this.RHO_GLOBAL* (this.maxAntEnergy / curBestEnergy + this.INIT_PHEROMONE / fitnessArry[i]);
//	//					System.out.print("vm: " + vm +  " pm: " + i);
//	//					System.out.println(" update pheromone:" + (pheromone[(int)vm][i] - bfPhr));
//					}
//				}
//			}
//			
//		}//end of if	
		
		for(int i = 0; i < numVM; i++)
		{
			int pm = placement[i];
			
//			if(!updatedPM[pm])
//			{
	//			double bfPhr = pheromone[i][pm];
				pheromone[i][pm] = (1 - this.RHO_GLOBAL) * pheromone[i][pm] + this.RHO_GLOBAL* (this.maxAntEnergy / curBestEnergy );
	//			System.out.println("update pheromone:" + (pheromone[i][pm] - bfPhr));
//			}
		}
		
	}// end of updateGlobalPheromone
	private double evaluateSinglePM(ArrayList server,int pm)
	{
		if(server.size() < 1)
			return Double.MAX_VALUE;
		
		double fitness = 0;
		double useCPU = 0;
		double useMEM = 0;
		double traffic = 0;
		double part1;
		double part2;
		int count = 0;
		int size = server.size();
		boolean [][]calTraffic = new boolean[size][size];
		
		for(int i = 0; i < size; i++)
			for(int j = 0; j < size; j++)
				calTraffic[i][j] = false;
		
		for(Object vm: server)
		{
			useCPU += this.vmCPU[(int)vm];
			useMEM += this.vmMEM[(int)vm];
		}
		
		part1 = Math.abs((this.pmCPU[pm] - useCPU) - (this.pmMEM[pm] - useMEM));
		part2 = (useCPU + useMEM);
		
		for(int i = 0 ; i < size; i++)
			for(int j =0; j < size; j++)
			{
				int vm1 = (int) server.get(i);
				int vm2 = (int) server.get(j);
				
				if(vm1 != vm2 && !calTraffic[i][j])
				{
					traffic += this.c[(int)vm1][(int)vm2];
					calTraffic[i][j] = true;
					calTraffic[j][i] = true;
					count ++;
				}
			}
		fitness = (part1 / part2 );
		
		if(traffic != 0)
		{
			fitness =  fitness + (1 - (traffic) / (count * this.maxTraffic)) / 10;
//			System.out.println();
//			System.out.println("maxTraffic:"+ this.maxTraffic * count);
//			System.out.println("traffic:" + traffic);
//			System.out.println("traffic:" +  (1 - (traffic) / (count * this.maxTraffic)));
//			System.out.println("part/part1 " + (part1 / part2 ) * ((1 - (traffic) / (count * this.maxTraffic)) / 10));
//			System.out.println();
//			
		}
//			
		return fitness;	
	}

	private boolean isResouceSatisfy(int vm, int pm, int[] useCPU, int[] useMEM)
	{
		if( useCPU[pm] + this.vmCPU[vm] < this.pmCPU[pm] &&
			useMEM[pm] + this.vmMEM[vm] < this.pmMEM[pm] )
			return true;
		
//		System.out.println("aco resource is not enough");
		return false;
	}// end of isResouceSatisfy
	private boolean isNetworkSatisfy(int index, int pm1, int[] placement)
	{
		int vm1 = this.vmList.get(index);
		
		for(int i = 0 ; i < index; i++)
		{
			int vm2 = this.vmList.get(i);
			int pm2 = placement[vm2];
			
			if(this.isLinkOverLoad(pm1, pm2, this.c[vm1][vm2]))
			{
//				System.out.println("aco network is not enough");
				return false;
			}
			
		}
		
		return true;
	}// end of isNetworkSatisfy
	
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
		
//			
//		if(this.isPELinkOverLoad(edge1, pod1, pm1, traffic) || this.isPELinkOverLoad(edge2, pod2, pm2, traffic))
//			return true;
		
		if(pod1 == pod2 && edge1 == edge2)
			return false;
		
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
	
	private void updateHeuristic(ArrayList serverList,int index, double[][] heuristic, int[] useCPU, int[] useMEM, boolean[] usePM, int[] placement) 
	{
		int vm = this.vmList.get(index);
		double tPr = this.BIG_FLOW_NUM / ((index + 1)) ;
	
		for(Object server:serverList)
		{
			int pm = (int)server;
			
			double pr = (index + 1 ) /(this.numVM * 1.0);
			double heuristicOne =  this.heuristicOne(pm, vm, useCPU, useMEM, usePM) ;// resource balance.
			double heuristicTwo =  this.heuristicTwo(pm, vm, useCPU, useMEM, usePM, placement);// server energy.
			
//			double heuristicFour = this.heuristicFour(placement);// switch energy.
//			double heuristicFive = this.heuristicFive(pm,index,placement);// traffic.
//			double heuristicSix  = this.heuristicSix(pm);
			
			
			
			heuristic[vm][pm] = ((1 + pr) * heuristicOne + (1 - pr) * heuristicTwo); // + heuristicFour) ;
//			heuristic[vm][pm] = heuristic[vm][pm] * this.heuristicThree(pm, vm, heuristic[vm][pm]);
			
	
//			System.out.println();
		}

	}//end of updateHeuristic
	
	private double heuristicOne(int pm,int vm, int[] useCPU, int[] useMEM,boolean[] usePM)// resource balance.
	{
		double sum = 0;
		double rCPU,rMEM;
		
		for(int i = 0; i < numPM; i++)
		{
			if(i == pm)
			{
				rCPU = this.pmCPU[i] - useCPU[i] - this.vmCPU[vm];
				rMEM = this.pmMEM[i] - useMEM[i] - this.vmMEM[vm];
				sum += Math.abs(rCPU - rMEM) / (useMEM[i] + this.vmMEM[vm] + useCPU[i] + this.vmCPU[vm] );	
			}
			else
			{
				rCPU = this.pmCPU[i] - useCPU[i] ;
				rMEM = this.pmMEM[i] - useMEM[i] ;
				
				if(useMEM[i] !=0 || useCPU[i] != 0 )
					sum += Math.abs(rCPU - rMEM) / (useMEM[i]  + useCPU[i]);	
			}
		}

		return 1/sum;
	}// end of heuristicOne
	private double heuristicTwo(int pm,int vm, int[] useCPU, int[] useMEM,boolean[] usePM, int[] placement)// server energy.
	{
		int[] pmWorkCPU = new int[this.numPM];
		boolean [] pmEmpty = new boolean[numPM];
		double total = 0,u;
		int count = 0;
		
		for(int i = 0; i < numPM; i++)
		{
			pmWorkCPU[i] = 0;
			pmEmpty[i] = true;
		}
		
		for(int i = 0; i < numVM; i++)
		{
			if(placement[i] != -1)
			{
				pmEmpty[placement[i]] = false;
				count ++;
			}
		}
		
		//calculate the total energy of pms.
		for(int i = 0; i < placement.length; i++)
			if(placement[i] != -1)
				pmWorkCPU[placement[i]] += this.vmCPU[i]; 
				
		pmWorkCPU[pm] += this.vmCPU[vm];
		
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
		
		if(total == 0)
			return 0;
		else
		{
//			System.out.println("heuristicTwo:" + (1/(total/count)));
			return (1/(total/count));
		}
		
		
	}// end of heuristicTwo
	
	private double heuristicThree(int pm,int vm, double exp) // traffic
	{
		
		double maxTraffic = Double.MIN_VALUE;
		int maxIndex = -1;
		
		for(int i = 0; i < vm ; i++)
		{
			if(maxTraffic < this.c[vm][i])
			{
				maxTraffic = this.c[vm][i];
				maxIndex = i;
			}
		}
		
//		System.out.println(heuristicOne * 100);

		if(pm == maxIndex)
		{
			return( 1 + maxTraffic / this.maxTraffic );
		}
		
		return 1;
	}// end of heuristicThree
	
	private double heuristicFour(int[] placement)// switch energy.
	{
	
		double total = 1 / this.totalSwitchEnergy(placement) * 100;
		
//		System.out.println("heuristicFour:" + total);
		return total ;

	}// end of heuristicFour
	
	private double heuristicFive(int pm,int index,int[] placement)// traffic .
	{
	
		int throughSwitch;
		double sumTraffic = 0;

		int vm1 = this.vmList.get(index);
		
		for(int i = 0; i < index ; i++)
		{		
			int vm2 = this.vmList.get(i);
			throughSwitch = this.throughSwitch[pm][placement[vm2]];
			sumTraffic += throughSwitch * this.c[vm1][vm2];
		}
		
		if(sumTraffic > 0)
			return 1 / sumTraffic ;	
		
		return 0;

	}// end of heuristicFive
	
	private double heuristicSix(int pm)// traffic .
	{
	
		int podNum  = this.getPodNum(pm);
		int edge = this.getEdge(podNum, pm);
		int pmIndex = this.getPmIndex(podNum, edge, pm);
		
		double linkUtilization = this.pod[podNum].pLinkE[pmIndex][edge] / this.PM_EDGE_CAPACITY;
		
		return Math.exp((Math.pow(linkUtilization, 10) / 0.04));
	

	}// end of heuristicFive
	
	
	public double getSumTraffic(int[] placement)
	{
		int throughSwitch;
		double sumTraffic = 0;
		
		for(int i = 0; i < numVM; i++)
			for(int j = 0; j < numVM; j++)
			{
//				System.out.println(j + " " + placement[j]);
				if(placement[i] != -1 && placement[j] != -1)
				{
					throughSwitch = this.throughSwitch[placement[i]][placement[j]];
					sumTraffic += throughSwitch * this.c[i][j];
				}
			}
		
		
//		System.out.println("aco totalTraffic:" + vmTotalTraffic);
//		System.out.println("aco:" + sumTraffic);
		
		return sumTraffic ;
			
	}// end of  bandWidth
	
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
		boolean useNewLink = true;
		
		boolean p = false;
		
		Pod pod = this.pod[podNum];
		
		for(int i = 0; i < pod.eLinkA[edge].length;i++)
		{
			if(pod.eLinkA[edge][i] < minTraffic)// &&  pod.eLinkA[edge][i] > 0)
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
			if(pod.aLinkC[aggre][i] < minTraffic)// && pod.aLinkC[aggre][i] > 0)
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
		this.coreUsed[aggre * (this.fatTreePod / 2) + core] = true;
		
		return core;
	}//end of getCore
	private int usePodNum(int pm)
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

		this.pod[pod].used = true;
		return pod;
	}//end of usePodNum
	
	private double weightArraySum(double [] weightArrays) {  
		double weightSum = 0;  
		for (double weightValue : weightArrays) {  
			weightSum += weightValue;  
		}  
		return weightSum;   
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
	
	private double getMeanMaxServerEnergy(int [] placement)
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
		
		int count = 0;
		for(int i = 0; i < numPM; i++)
		{
			if(!pmEmpty[i])
			{
				count++;
				
				if(pmCPU[i] > this.SERVER2_CPU_MIPS)
					total += this.HIGH_SERVER_FULL_ENERGY;
				else
					total += this.SERVER_FULL_ENERGY;
			}
		}
			
		return total / count;
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
	

	public double getAcoTotalEnergy()
	{
		return this.totalEnergy;
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
	
	public static void generateSence(int vm,int pm) throws IOException
	{
		int VM = vm,PM = pm;
		int testCount = 0;
	
		int ffd_v = 0;
		int ff_v = 0;
		int bfd_v = 0;
		int bf_v = 0;
		int aco_v =0;
		int ga_v = 0;
	

		
		String fileName = "D:\\requirement\\" + VM + "-" + PM + "\\";
		HashMap<String , Object> parameters = new HashMap<String , Object>();
		parameters.put("ants",20);
		parameters.put("maxGen",10);
		parameters.put("q",0.95);
		
		int []useCPU = new int[PM];
		int []useMEM = new int[PM];
		
		List<String> fileList = GetFoldFileNames.getFileName(fileName + "TRAFFIC");
		
		int runTimes = 0;
		int betterTimes = 0;
		int  fesibleRunTimes = 0;
		
		String file = fileList.get(0);
		String name = file;
		for(;; ) {
			double saveEnergy = 0;
			double saveTraffic = 0 ;
			
			runTimes ++;
			int [] aco_placement;
			int [] ffd_placement;
			int [] ff_placement;
			int [] acosd_placement = null;
			int [] bfd_placement;
			int [] bf_placement;
			
			double aco_bw = -1;
			double ga_bw = -1;
		
			FfdEVMP ffd_evmp = new FfdEVMP(VM,PM);
			FfEVMP  ff_evmp = new  FfEVMP(VM,PM);

			AcoEVMP aco_evmp = new AcoEVMP(VM,PM,fileName,file);
//			AcoSdEVMP  acosd_evmp = new AcoSdEVMP(VM,PM);

//			acosd_evmp.getResource(aco_evmp.getPmCPU(), aco_evmp.getPmMEM(), aco_evmp.getVmCPU(), aco_evmp.getVmMEM(), aco_evmp.getVmTraffic());
			ffd_evmp.getResource(aco_evmp.getPmCPU(), aco_evmp.getPmMEM(), aco_evmp.getVmCPU(), aco_evmp.getVmMEM(), aco_evmp.getVmTraffic());
			ff_evmp.getResource(aco_evmp.getPmCPU(), aco_evmp.getPmMEM(), aco_evmp.getVmCPU(), aco_evmp.getVmMEM(), aco_evmp.getVmTraffic());
			System.out.println("run:" + runTimes + " ----start simulation #VM-" + VM + " #PM-" + PM  + " " + file + "-----");
			
	
			ffd_placement = ffd_evmp.FFD();
			ff_placement = ff_evmp.FF();

//			acosd_placement = acosd_evmp.ACOSD(parameters);
			
			long startTime = System.currentTimeMillis(); 
			aco_placement = aco_evmp.ACO(parameters);
			long endTime = System.currentTimeMillis(); 
			double acoTime = (endTime - startTime);
			
			
//			aco_evmp.placeDetil("aco:", aco_placement);
//			aco_evmp.placeDetil("ffd:", ffd_placement);
			
			
			java.text.DecimalFormat df = new java.text.DecimalFormat("#.00");  
			
			System.out.println("ACO-Time:"+ acoTime/1000.0 + " s");
	
			if(aco_placement == null)
			{
				continue;
			}
			else {
				saveEnergy += aco_evmp.getAcoTotalEnergy();
				System.out.println("ACO-TotalEnergy:"+Double.valueOf(df.format(aco_evmp.getAcoTotalEnergy())));
				aco_bw = Double.valueOf(df.format(aco_evmp.getSumTraffic(aco_placement)));
				
				saveTraffic += aco_bw;
				
				System.out.println("ACO :" + aco_bw );
				System.out.println();
			}
			
//			if(acosd_placement == null)
//			{
////				ga_v ++;	
//			}
//			else {
//				
////				saveEnergy = acosd_evmp.getAcoSdTotalEnergy() - saveEnergy;
//				
//				System.out.println("acosd-TotalEnergy:"+Double.valueOf(df.format(acosd_evmp.getAcoSdTotalEnergy())));
//				ga_bw = Double.valueOf(df.format(aco_evmp.getSumTraffic(acosd_placement)));
//				saveTraffic =  ga_bw - saveTraffic;
//				System.out.println("aco-sd:" + ga_bw );
//				System.out.println();
//			}
			
			if(ffd_placement == null)
			{
				ffd_v ++;
			}
			else
			{
				System.out.println("FFD-TotalEnergy:"+Double.valueOf(df.format(ffd_evmp.getFfdTotalEnergy())));
				ffd_evmp.ffd_bw = Double.valueOf(df.format(aco_evmp.getSumTraffic(ffd_placement)));
				System.out.println("FFD-useMeanSwitch:" + ffd_evmp.ffd_bw );
				System.out.println();
			}
			
			if(ff_placement == null)
			{
				ff_v ++;
			}
			else
			{
//				System.out.println("FF-TotalEnergy:"+Double.valueOf(df.format(ff_evmp.getFfTotalEnergy())));
				ff_evmp.ff_bw = Double.valueOf(df.format(aco_evmp.getSumTraffic(ff_placement)));
//				System.out.println("FF-useMeanSwitch:" + ff_evmp.ff_bw + " switches");
//				System.out.println();
			}
			
			
			
			System.out.println("Demand BandWidth:" + df.format(aco_evmp.vmTotalTraffic) + " MB");
			System.out.println();
		
			if(aco_placement != null && ffd_placement != null )
			{
				fesibleRunTimes++;
				System.out.println("save traffic :" + (ffd_evmp.ffd_bw - aco_bw) );
		
				if((ffd_evmp.getFfdTotalEnergy() - aco_evmp.getAcoTotalEnergy() ) > 0 && (ffd_evmp.ffd_bw - aco_bw) > 0) 
					betterTimes ++;
			}
	
			if(aco_placement != null && ffd_placement != null )
				System.out.println("save energy:" +(ffd_evmp.getFfdTotalEnergy() - aco_evmp.getAcoTotalEnergy() ) );
		
			System.out.println();
//			System.out.println("FFD: violation " + ffd_v / (runTimes * 1.0) );
//			System.out.println("FF: violation " +  ff_v / (runTimes * 1.0) );
//			System.out.println("ACO: violation " + aco_v / (runTimes * 1.0) );
//			System.out.println("GA: violation " + ga_v / (runTimes * 1.0) );
			System.out.println("fesible run:" + fesibleRunTimes  +"  better:" + betterTimes);
	
			if(     ffd_evmp.getFfdTotalEnergy() - aco_evmp.getAcoTotalEnergy() > 0 && 
					ff_evmp.getFfTotalEnergy() - aco_evmp.getAcoTotalEnergy() > 0 && 
					
//					ga_evmp.getGaTotalEnergy() - aco_evmp.getAcoTotalEnergy() > 0 &&
					
					ffd_evmp.ffd_bw - aco_bw > 0.5 &&
					ff_evmp.ff_bw - aco_bw > 0 )
			{
				String time = new SimpleDateFormat("yyyyMMddHHmmss").format(new Date());
			
				RequirementWR.generateRequirement(aco_evmp.getPmCPU(),PM, 10,  "D:\\REQUIREMRNT\\"+ VM + "-" + PM + "\\PM-CPU\\"    + (ffd_evmp.getFfdTotalEnergy() - aco_evmp.getAcoTotalEnergy()) + name);
				RequirementWR.generateRequirement(aco_evmp.getPmMEM(),PM, 10, "D:\\REQUIREMRNT\\"+ VM + "-" + PM + "\\PM-MEM\\"  + (ffd_evmp.getFfdTotalEnergy() - aco_evmp.getAcoTotalEnergy()) + name);
 			    RequirementWR.generateRequirement(aco_evmp.getVmMEM(),VM, 10,  "D:\\REQUIREMRNT\\"+ VM + "-" + PM + "\\VM-CPU\\"   + (ffd_evmp.getFfdTotalEnergy() - aco_evmp.getAcoTotalEnergy()) + name);
				RequirementWR.generateRequirement(aco_evmp.getVmCPU(),VM, 10,"D:\\REQUIREMRNT\\"+ VM + "-" + PM + "\\VM-MEM\\"    + (ffd_evmp.getFfdTotalEnergy() - aco_evmp.getAcoTotalEnergy()) + name);
				RequirementWR.generateRequirement(aco_evmp.getTraffic(),VM, VM, "D:\\REQUIREMRNT\\"+ VM + "-" + PM +"\\TRAFFIC\\"   + (ffd_evmp.getFfdTotalEnergy() - aco_evmp.getAcoTotalEnergy()) + name);
				
			}
			
			System.out.println();
		}
		
	}
	public static void testAco(int vm,int pm) throws IOException
	{
		int VM = vm,PM = pm;
		int testCount = 30;
	
		int ffd_v = 0;
		int ff_v = 0;
		int bfd_v = 0;
		int bf_v = 0;
		int aco_v =0;
		int ga_v = 0;
	
		
		String fileName = "D:\\requirement-new\\" + VM + "-" + PM + "\\";
		HashMap<String , Object> parameters = new HashMap<String , Object>();
		parameters.put("ants",20);
		parameters.put("maxGen",10);
		parameters.put("q",0.99);
		
		int []useCPU = new int[PM];
		int []useMEM = new int[PM];
		
		List<String> fileList = GetFoldFileNames.getFileName(fileName + "TRAFFIC");
		
		int [] placement;
		int runTimes = 0;
		int betterTimes = 0;
		int  fesibleRunTimes = 0;
		java.text.DecimalFormat df = new java.text.DecimalFormat("#.00");  
		
		double [] energy = new double [testCount];
		double [] useTraffic = new double[testCount];
		double [] runTime = new double[testCount];
		double[] switchchassis = new double[1];
		
		
		
		for(String file : fileList ) {
			
			System.out.println("Start " + VM + "-" + PM + file + " test");

			int i = 0;
			while(i < testCount)
			{
				placement = null;
				double aco_bw = -1;
				
				AcoEVMP aco_evmp = new AcoEVMP(VM,PM,fileName,file);
//				PSOBUPT evmp = new PSOBUPT(VM,PM,200);
//				evmp.getResource(aco_evmp.getPmCPU(), aco_evmp.getPmMEM(), aco_evmp.getVmCPU(), aco_evmp.getVmMEM(), aco_evmp.getTraffic());
				
				switchchassis[0] = aco_evmp.switchChassis;
				System.out.println("vm:" + VM + " pm:" + PM + " :" + aco_evmp.fatTreePod );
//				i++;
				
				long startTime = System.currentTimeMillis(); 
				placement = aco_evmp.ACO(parameters);
				long endTime = System.currentTimeMillis(); 
				long acoTime = (endTime - startTime);
				

				if(placement != null)
				{
					energy[i] = Double.valueOf(df.format(aco_evmp.getAcoTotalEnergy()+ aco_evmp.switchChassis));
					useTraffic[i] = Double.valueOf(df.format(aco_evmp.getSumTraffic(placement)));
					runTime[i] = acoTime/1000.0;
					System.out.println(i);
					System.out.println("energy:" + energy[i]);
					System.out.println(file + " vm:" + VM + " pm:" + PM + " :" + "runTime:" + runTime[i]);
					System.out.println();
					i++;	
				}
			}
			
			double[] average = new double[1];
			double[] StandardDevition = new double[1];
			
			average[0] = TestResultWR.getAverage(energy, testCount);
			StandardDevition[0] =  TestResultWR.getStandardDevition(energy, testCount, average[0]);
			
			int ants = (int) parameters.get("ants");
			int gen = (int) parameters.get("maxGen");
//			
			String al = "NEW-ACO";
			
			TestResultWR.generateResult("energy", energy, testCount, 5, "D:\\scenesNew\\" + al + "\\" + file);
			TestResultWR.generateResult("useTraffic", useTraffic, testCount, 5, "D:\\scenesNew\\"+ al + "\\" + file);
			TestResultWR.generateResult("runTime", runTime, testCount, 5, "D:\\scenesNew\\" + al + "\\" + file);
			
			TestResultWR.generateResult("energy mean", average, 1, 5, "D:\\scenesNew\\" + al + "\\" + file);
			TestResultWR.generateResult("energy standardDevition", StandardDevition, 1, 5, "D:\\scenesNew\\"+ al + "\\" + file);
//			
			average[0] = TestResultWR.getAverage(useTraffic, testCount);
			StandardDevition[0] =  TestResultWR.getStandardDevition(useTraffic, testCount, average[0]);
			
			TestResultWR.generateResult("useTraffic mean", average, 1, 5, "D:\\scenesNew\\"+ al + "\\" + file);
			TestResultWR.generateResult("useTraffic standardDevition", StandardDevition, 1, 5,"D:\\scenesNew\\"+ al + "\\" + file);
			
			average[0] = TestResultWR.getAverage(runTime, testCount);
			
			TestResultWR.generateResult("runTime mean", average, 1, 5, "D:\\scenesNew\\"+ al + "\\" + file);
//			TestResultWR.generateResult("useTraffic standardDevition", StandardDevition, 1, 5,"D:\\scenesNew\\"+ al + "\\" + file);
			
			TestResultWR.generateResult("switchChassis", switchchassis, 1, 5,"D:\\scenesNew\\"+ al + "\\" + file);
			
			
			System.out.println();
		}
		
	}
	
	
	public static void generateDataSet(int vm,int pm) throws IOException
	{
		int VM = vm,PM = pm;
	
		String fileName = "D:\\requirement-new\\" + VM + "-" + PM + "\\";
		HashMap<String , Object> parameters = new HashMap<String , Object>();
		parameters.put("ants",20);
		parameters.put("maxGen",10);
		parameters.put("q",0.99);
		
		List<String> fileList = GetFoldFileNames.getFileName(fileName + "TRAFFIC");
		
		
		for(String file : fileList ) {
			
			System.out.println("Start " + VM + "-" + PM + file + " test");
			AcoEVMP aco_evmp = new AcoEVMP(VM,PM,fileName,file);
			
			String N_Eflow_Cluster = null;
			String N_Mflow_Cluster = null;
			
			int minEflow = 0 ;
			int maxEflow = 0 ;
			
			switch(pm)
			{
			case 20:
				N_Eflow_Cluster = "[5,10]";
				
				minEflow = 5;
				maxEflow = 10;
				
				N_Mflow_Cluster = "[3,10]";
				break;
			case 40:
				N_Eflow_Cluster = "[5,15]";
				minEflow = 5;
				maxEflow = 15;
				N_Mflow_Cluster = "[5,10]";
				break;
			case 60:
				N_Eflow_Cluster = "[10,15]";
				minEflow = 10;
				maxEflow = 15;
				N_Mflow_Cluster = "[5,15]";
				break;
			case 80:
				N_Eflow_Cluster = "[15,20]";
				minEflow = 15;
				maxEflow = 20;
				N_Mflow_Cluster = "[10,15]";
				break;
			
			}
			
			
			String fn = "D:\\DataSet\\" + file + "\\" ;
			
			
			//file -- Summarize.txt
			File f = new File(fn + "1. Introduction.txt");  
			FileWriter out = new FileWriter(f,true);  
			
			out.write("Please note that: \r\n");
			out.write("* the cloud datacenter is based on the commonly used topology, Fat-Tree.\r\n");
			out.write("* There are 2 types of PMs. Each PM is associated with CPU and MEM capacities,and POWER. See file Configuration_PMs.txt for details.\r\n");
			out.write("* There are 2 types of switches, one for core level switching and the other for aggregation/edge level switching. Each Switch is associated "
					+ "with POWER_Chassis and POWER_Port. POWER_Chassis is the power consumption from the chassis and POWER_Port is the power consumption from a switch port."
					+ " Note that, a port with different maximum data rates consumes different amounts of power. So, POWER_Port is associated with maximum data rate (MAX_DATA_RATE)."
					+ " See file Configuration_Switches.txt for details.\r\n");
			
			out.write("* Each VM is associated with CPU and MEM requirements. See file PhysicalResourceRequirement_VMs.txt for details.\r\n");
			out.write("* Traffic demand between each VM pair is in file Traffic_VM_Pair.txt.\r\n");
			out.write("\r\n");
			out.write("\r\n");
			out.write("The following summarizes this instance:\r\n");
			out.write("\r\n");
			
			out.write("Number of VMs:\r\n");
			out.write("|V|=" + VM + "\r\n");
			out.write("\r\n");
			
			out.write("Number of PMs:\r\n");
			out.write("|P|=" + PM + "\r\n");
			out.write("\r\n");
			
			out.write("Number of core switches:\r\n");
			out.write("N_Switch_Core =" + aco_evmp.core + "\r\n");
			out.write("\r\n");
			
			out.write("Number of aggregation  switches:\r\n");
			out.write("N_Switch_Core =" + aco_evmp.aggre + "\r\n");
			out.write("\r\n");
			
			out.write("Number of edge  switches:\r\n");
			out.write("N_Switch_Core =" + aco_evmp.edge + "\r\n");
			out.write("\r\n");
			
			out.write("Cluster size of an Elephant flow (E-flow) is randomly generated in the range:\r\n");
			out.write(N_Eflow_Cluster  + "\r\n");
			out.write("\r\n");
			
			out.write("Traffic demand between each VM pair in E-Flow is randomly generated in the range:\r\n");
			out.write("[25,100] Mb/s"  + "\r\n");
			out.write("\r\n");
			
			out.write("Cluster size of a Mouse flow (M-flow) is randomly generated in the range:\r\n");
			out.write(N_Mflow_Cluster  + "\r\n");
			out.write("\r\n");
			
			out.write("Traffic demand between each VM pair in M-Flow is randomly generated in the range:\r\n");
			
			if(VM / PM > 3)
				out.write("[1,3] Mb/s"  + "\r\n");
			else
				out.write("[1,5] Mb/s"  + "\r\n");
			out.write("\r\n");
			
			out.close();
			
//			System.out.println("|V|\t" + "|P|\t" + "N_Switch_Core\t" + "N_Switch_Aggregation\t" + "N_Switch_Edge\t" + "N_Eflow_Cluster\t" + "N_Mflow_Cluster\t");
//			System.out.println(VM +"\t" + PM+"\t    " + aco_evmp.core+"\t\t      " + aco_evmp.aggre+"\t\t    " + aco_evmp.edge+"\t\t    " + N_Eflow_Cluster+"    " +  N_Mflow_Cluster + "\t") ;
			
			
			//file -- vm config
			f = new File(fn + "4. PhysicalResourceRequirement_VMs.txt");  
			out = new FileWriter(f,true);
			out.write("* This file specifies the physical resource requirement for placing VMs. \r\n");
			out.write("\r\n");
//			System.out.println("VM\t" + "CPU\t" +"MEM\t");
			out.write("VM\t" + "CPU(MIPS)\t" +"MEM(MB)\t" + "\r\n");
			for(int i = 0 ; i < VM; i++ )
			{
//				System.out.println(i+"\t" + aco_evmp.getVmCPU()[i]+"\t" + aco_evmp.getVmMEM()[i]+ "\t");
				out.write(i+"\t" + aco_evmp.getVmCPU()[i]+"         \t" + aco_evmp.getVmMEM()[i]+ " \t" + "\r\n");
			}
			out.close();
			
			//file -- pm config
			f = new File(fn + "2. Configuration_PMs.txt");  
			out = new FileWriter(f,true);
			out.write("* This file specifies the configuration of PMs. \r\n");
			out.write("* There are two types of PMs in use, namely Type-1 and Type-2. \r\n");
			out.write("\r\n");
			
			
//			System.out.println("PM\t" + "CPU\t" +"MEM\t" + "ENERGY\t");
			out.write("PM\t" + "TYPE\t\t" + "CPU(MIPS)\t" +"MEM(MB)\t\t" + "Energy(W)\t" +"\r\n");
			for(int i = 0 ; i < PM; i++ )
			{
				if(aco_evmp.getPmCPU()[i] > SERVER2_CPU_MIPS )
					out.write(i+"\t" +"Type-1\t\t" +aco_evmp.getPmCPU()[i]+"         \t" + aco_evmp.getPmMEM()[i]+ "         \t" + HIGH_SERVER_FULL_ENERGY + " \t" +"\r\n");
				else
					out.write(i+"\t" + "Type-2\t\t" +aco_evmp.getPmCPU()[i]+"         \t" + aco_evmp.getPmMEM()[i]+ "         \t" + SERVER_FULL_ENERGY + " \t" + "\r\n");
					
			}
			
			out.close();
			
			//file -- traffic
			
			f = new File(fn + "5. Traffic_VM_Pair.txt");  
			out = new FileWriter(f,true);
			
			double [][] t = aco_evmp.getTraffic();
			
			out.write("* This file specifies the traffic demands among VMs. \r\n");
			out.write("* Two flows are considered, namely E-flow (elephant flow) and M-flow (mouse flow). Values larger than 25 Mb/s are regarded as E-flows. \r\n");
			out.write("\r\n");
			
			out.write("VM_i\t" + "VM_j\t" +"Traffic(Mb/s)\t" + "\r\n");
			
			int nEflow = 0;
			for(int i = 0; i < VM; i++)
				for(int j = 0; j < VM; j++)
				{
					if(t[i][j] > 0)
					{
						if(t[i][j] > 25)
							nEflow ++;
						
						while(t[i][j] > 100)
							t[i][j] -= 1;
						
						if(VM / PM > 4)
							while(t[i][j] < 25 && t[i][j] >3)
								t[i][j] -= 1;
						else
							while(t[i][j] < 25 && t[i][j] >5)
								t[i][j] -= 0.5;
					}
				}
			
			int num = minEflow + new Random().nextInt(maxEflow - minEflow);
			nEflow = nEflow / 2;
			
			while(num - nEflow > 0)
			{
				
				int  i = new Random().nextInt(VM);
				
				for(int j = 0; j < VM; j++)
					if(t[i][j] == 0)
					{
						t[i][j] = 25 + new Random().nextInt(100 - 25);
						t[j][i] = t[i][j];
						nEflow ++;
						
						break;
					}
				
			}
			
			java.text.DecimalFormat df = new java.text.DecimalFormat("0.00");  
			
			for(int i = 0; i < VM; i++)
				for(int j = 0; j < VM; j++)
				{
					if(t[i][j] > 0)
					{
						out.write(i+"\t" + j+"\t" +df.format(t[i][j]) + " \t" + "\r\n");
					}
				}
					
			out.close();
			
			//file -- switches
			f = new File(fn + "3. Configuration_Switches.txt");  
			out = new FileWriter(f,true);
			out.write("* This file specifies the configuration of Switches. \r\n");
			out.write("* There are two types of switches in use. type_core is for core level switching and type_agg,edge is for aggregation and edge level switching.\r\n");
			out.write("\r\n");
			
			//Core Switch 10MB 100M 1000MB
			out.write("LEVEL\t"+ "TYPE\t\t" + "POWER_Chassis(W)\t" + "MAX_DATA_RATE(Mb/s)\t" + "POWER_Port(W)\t" + "\r\n");
			for(int i = 0; i <1; i++)
			{
				out.write("Core\t" +"type_core\t  "+ CORE_SWITCH_IDLE_ENERGY+"     \t\t " + "10   \t\t\t" + _10MB_ENERGY * 20 + " \t   " + "\r\n");
				out.write("Core\t" +"type_core\t  "+ CORE_SWITCH_IDLE_ENERGY+"     \t\t " + "100   \t\t\t" + _100MB_ENERGY *20 +" \t   " + "\r\n");
				out.write("Core\t" +"type_core\t  "+ CORE_SWITCH_IDLE_ENERGY+"     \t\t "	 + "1000  \t\t\t" + _1000MB_ENERGY *20 +" \t   " + "\r\n");
			}
			
			for(int i = 0; i < 1; i++)
			{
				out.write("Agg\t" +"type_agg\t  "+ EDGE_SWITCH_IDLE_ENERGY+"     \t\t " + "10   \t\t\t" + _10MB_ENERGY  + " \t   " + "\r\n");
				out.write("Agg\t" +"type_agg\t  "+ EDGE_SWITCH_IDLE_ENERGY+"     \t\t " + "100   \t\t\t" + _100MB_ENERGY  +" \t   " + "\r\n");
				out.write("Agg\t" +"type_agg\t  "+ EDGE_SWITCH_IDLE_ENERGY+"     \t\t "	 + "1000  \t\t\t" + _1000MB_ENERGY +" \t   " + "\r\n");
			}
			
			for(int i = 0; i < 1; i++)
			{
				out.write("Edge\t" +"type_edge\t  "+ EDGE_SWITCH_IDLE_ENERGY+"     \t\t " + "10   \t\t\t" + _10MB_ENERGY  + " \t   " + "\r\n");
				out.write("Edge\t" +"type_edge\t  "+ EDGE_SWITCH_IDLE_ENERGY+"     \t\t " + "100   \t\t\t" + _100MB_ENERGY  +" \t   " + "\r\n");
				out.write("Edge\t" +"type_edge\t  "+ EDGE_SWITCH_IDLE_ENERGY+"     \t\t "	 + "1000  \t\t\t" + _1000MB_ENERGY +" \t   " + "\r\n");
			}
				
			out.close();
			
			System.out.println();
			
		
		}
		
	}
	
	
	
	
	
	public static void main(String[] args) throws IOException
	{
		
		generateDataSet(40,20);
		generateDataSet(60,20);
		generateDataSet(120,20);
		
		generateDataSet(80,40);
		generateDataSet(120,40);
		generateDataSet(240,40);
		
		generateDataSet(120,60);
		generateDataSet(180,60);
		generateDataSet(360,60);
		
		generateDataSet(160,80);
		generateDataSet(240,80);
		generateDataSet(480,80);
		

	}// end of main.	
	
	
	
	
	
}
