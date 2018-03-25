package VMP_Algorithm;

import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

public class AcoSdEVMP {
	

	public double acosd_bw = 0;
	
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

	private final static double SERVER2_CPU_MIPS = 24000;
	private final static double SERVER2_MEM = 16000;
	
	private final static double HIGH_SERVER_FULL_ENERGY = 550;//252;
	private final static double SERVER_FULL_ENERGY = 420;//252;
	private final static double SERVER_IDLE_ENERGY = SERVER_FULL_ENERGY * 0.7;
	private final static double HIGH_SERVER_IDLE_ENERGY = HIGH_SERVER_FULL_ENERGY * 0.7;
	
	
	private final static double SWITCH_IDLE_ENERGY = 147 ;
	
	private final static double CORE_SWITCH_FULL_ENERGY = 650 ;
	private final static double CORE_SWITCH_IDLE_ENERGY = 550 ;
	
	
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

	private final static double SERVER_CPU = 18080 * 2.4;
	private final static double SERVER_MEM = 64 * 1024;
	
	private final static double INIT_PHEROMONE = 0.0001;
	private final static double INIT_HEURISTIC = 0;
	private final static double GUID_HEURISTIC = 0.1;
	
	private final static double BETA = 3 ;
	private final static double ALPHA = 1;
	private final static double EPSILON = 0.0001;
	private final static double RHO_LOCAL = 0.35;
	private final static double RHO_GLOBAL = 0.35;
	private final static double VM_PROPORTION = 0.3; 
	
	private final static int    BIG_FLOW_NUM = 10;
	
	private final double PR_BIAS = 0.1; 


	public AcoSdEVMP(int numVm, int numPm) 
	{
		this.numPM = numPm;
		this.numVM = numVm;
		this.bestPlacement = new int[numVM];
		
		this.generateTopology(numPM);	
		this.generateThroughSwitch();

		
		for(int i = 0; i < numVM;i++)
		{
			
				this.vmList.add(i);
			
		}	
		
	}// end of ACO_EVMP
	
	
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
		 
		 this.generateVmList();
	}
	
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
		
		for(int i = 0; i < numVM; i++)
			vmList.add(i);
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
	public int [] ACOSD(HashMap<String , Object> parameters )
	{	
		int ants = 20;
		int maxGen = 10;
		double q = 0.95;
		
		
		double[] perAntEnergy = new double[ants];
		int [][]placeSet = new int [ants][numVM];
		double [][] pheromone = new double [numVM][numPM];
		double [][] heuristic = new double [numVM][numPM];
		
		int[] useCPU = new int[numPM];
		int[] useMEM = new int[numPM];
		boolean[] usePM  = new boolean[numPM];
		int[] bestPlacement;
		this.initTopology();
		this.initACO(pheromone,heuristic); // init parameter.
		int curGen = 0;
		int pm;
		
		while(curGen < maxGen)
		{
			for(int k = 0; k < ants; k++ )
			{
				// init the server's resoures used by previous ant.
				this.initServerUsedResource(useCPU,useMEM,usePM,placeSet[k],perAntEnergy);
				
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
			} // end of ants
			
			bestPlacement = this.getBestPlacement(placeSet,perAntEnergy);
			this.updateGlobalPheromone(bestPlacement,pheromone,curBestEnergy);
	
			curGen ++;
		} // end of while
		
		return this.bestPlacement;
	}// end of ACO
	
	
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
	
	private void calculateTotalEnergy(int[] placement, double[] perAntEnergy, int k) 
	{
		double pmEnergy = this.totalPmEnergy(placement);
		double switchEnergy = this.totalSwitchEnergy(placement);
		
	
//		System.out.println("ACO"  +  pmEnergy);
//		System.out.println("ACO"  +  switchEnergy);
		
		perAntEnergy[k] =  pmEnergy + switchEnergy;
		
		if(perAntEnergy[k] == 0)
			perAntEnergy[k] = Double.MAX_VALUE;
		
	}//end of calculateTotalEnergy
	
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
	private void initServerUsedResource(int[] useCPU, int[] useMEM, boolean[] usePM, int[] placement, double[] perAntEnergy)
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
		
		for(double energy : perAntEnergy) 
			energy = 0;
		
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
		double []exp = new double[serverList.size()];
		List<Integer> candidateServerIndex = new ArrayList<Integer>();
	
		
		for(int i = 0; i < serverList.size(); i++ )
		{	
			pm = (int) serverList.get(i);		

			exp[i] = pheromone[vm][pm] * Math.pow(heuristic[vm][pm],BETA);

			if(maxExp < exp[i])
			{
				placement[vm] = (int) serverList.get(i);	
				maxExp = exp[i];
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
			pr[i] = pheromone[vm][pm] *  Math.pow(heuristic[vm][pm] ,BETA);	
		}
		
		WeightRandom random = new WeightRandom();
		double weightSum = weightArraySum(pr);
		
		int pm = (int) serverList.get(random.getWeightRandom(pr, weightSum));
		placement[vm] = pm;
		usePM[pm] = true;

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
		int switches = core + edge + aggre;
		
		for(int k = 0; k < placeSet.length; k++)
		{  
			int [] placement = placeSet[k];
			curEnergy = perAntEnergy[k];
			curUseSwitch = this.useMeanSwitch(placement);
			if(minEnergy > curEnergy && (curUseSwitch  - minUseSwitch) / minUseSwitch < 0.1)
			{
				minEnergy = curEnergy;
				minUseSwitch = curUseSwitch;
				System.arraycopy(placement, 0, bestPlacement, 0, numVM);
			}
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

		for(int i = 0; i < numVM; i++)
		{
			int pm = placement[i];
			double bfPhr = pheromone[i][pm];
			pheromone[i][pm] = (1 - this.RHO_GLOBAL) * pheromone[i][pm] + this.RHO_GLOBAL/ curBestEnergy ;
		}
		
	}// end of updateGlobalPheromone

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
			double heuristicOne =  this.heuristicOne(pm, vm, useCPU, useMEM, usePM) ;// resource balance.
			double heuristicTwo =  this.heuristicTwo(pm, vm, useCPU, useMEM, usePM, placement);// server energy.
			heuristic[vm][pm] = (heuristicOne + heuristicTwo ); // + heuristicFour) ;
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
				total += (this.SERVER_IDLE_ENERGY + (1 - 0.7) * this.SERVER_FULL_ENERGY * u) / SERVER_FULL_ENERGY ;
			}
		}
		
		if(total == 0)
			return 0;
		else
		{
//			System.out.println("heuristicTwo:" + (1/(total/count)));
			return 1/(total);
		}
		
		
	}// end of heuristicTwo
	
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
		
		
//		System.out.println("aco totalTraffic:" + vmTotalTraffic);
//		System.out.println("aco:" + sumTraffic);
		
		return sumTraffic / this.vmTotalTraffic;
			
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
	
	

	public double getAcoSdTotalEnergy()
	{
		return this.totalEnergy;
	}

}// end of class ACO_EVMP
