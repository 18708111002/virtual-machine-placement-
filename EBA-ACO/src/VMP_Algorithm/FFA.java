package VMP_Algorithm;

import java.io.IOException;
import java.util.ArrayList;


public class FFA {

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
	private double ff_energy;
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
	
	public FFA(int vm,int pm)
	{
		this.numVM = vm;
		this.numPM = pm;
		
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
	
	
	public int [] FFA()
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
		
		boolean [] vmPlaced = new boolean[numVM];
		
		for(int i = 0; i < numVM; i++)
			vmPlaced[i] = false;
		
		for(int i = 0 ; i < numPM; i++)
		{
		
			for(int j = 0; j < numVM ; j++)
			{
				if( this.isResouceSatisfy(j,i,useCPU,useMEM) &&
						this.isNetworkSatisfy(j,i,ffPlacement,vmPlaced) )
				{
					if(!vmPlaced[j])
					{
						vmPlaced[j] = true;
						ffPlacement[j] = i;
						this.updateLinkLoadByFF(j, i, ffPlacement,vmPlaced);
						this.updatePmResource(j, i, useCPU, useMEM);
					}
				}

			}
		}
		
		for(int i = 0; i < numVM; i++)
		{
			if(vmPlaced[i] == false)
			{
				System.out.println("no placement " + i);
			}
		}
		
		double pmEnergy = this.totalPmEnergy(ffPlacement);
		double switchEnergy = this.totalSwitchEnergy(ffPlacement);
		
		this.ff_energy = pmEnergy + switchEnergy;
		java.text.DecimalFormat df = new java.text.DecimalFormat("#.00");
		this.ff_bw = Double.valueOf(df.format(this.useMeanSwitch(ffPlacement)));
		
		return ffPlacement;	
	}// end of FFD.
	
	private void updatePmResource(int vm, int pm, int[] useCPU, int[] useMEM)
	{
		useCPU[pm] += this.vmCPU[vm];
		useMEM[pm] += this.vmMEM[vm];
		
	}// end of updatePmResource
	

	
	private boolean isResouceSatisfy(int vm, int pm, int[] useCPU, int[] useMEM)
	{
		if( useCPU[pm] + this.vmCPU[vm] < this.pmCPU[pm] &&
			useMEM[pm] + this.vmMEM[vm] < this.pmMEM[pm] )
			return true;
		
		return false;
	}// end of isResouceSatisfy
	private boolean isNetworkSatisfy(int vm1, int pm1, int[] placement, boolean[] vmPlaced)
	{
		for(int i = 0 ; i < numVM; i++)
		{
			if(vmPlaced[i])
			{
				int vm2 = i;
				int pm2 = placement[vm2];
				
				if(this.isLinkOverLoad(pm1, pm2, this.c[vm1][vm2]))
				{
	//				System.out.println("network is not enough");
					return false;
				}
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
	
	private void updateLinkLoadByFF(int vm,int pm1,int[] placement, boolean[] vmPlaced)
	{
		
		for(int i = 0; i < numVM; i++)
		{
			if(!vmPlaced[i])
				continue;
			
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
	
	public double getFfTotalEnergy()
	{
		return this.ff_energy;
	}//end of getFfdTotalEnergy

}
