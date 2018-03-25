package VMP_Algorithm;
import java.util.Random;

public class Pod
{
	private int    limit  = 1;
	public  double [][]eLinkA;
	public  double [][]aLinkC;
	public  double [][]pLinkE;
	public  boolean used;
	public  boolean []edgeUsed;
	public  boolean []aggreUsed;
	
	private double generateBackTraffic(int limit)
	{
		return limit - Math.random();
	}
	
	public Pod(int k)
	{
		this.used = false;
		this.edgeUsed = new boolean[k/2];
		this.aggreUsed = new boolean[k/2];
		this.eLinkA = new double[k/2][];
		this.aLinkC = new double[k/2][];
		this.pLinkE = new double[k/2][];
		
		for(int i = 0; i < k / 2; i++)
		{
			this.eLinkA[i] = new double[k/2];
			this.aLinkC[i] = new double[k/2];
			this.pLinkE[i] = new double[k/2];
			
			this.edgeUsed[i] = false;
			this.aggreUsed[i] = false;
			
			
			for(int j = 0; j < k / 2; j++)
			{
				this.eLinkA[i][j] = 0;//this.generateBackTraffic(limit);
				this.aLinkC[i][j] = 0;//this.generateBackTraffic(limit);
				this.pLinkE[i][j] = 0;
			}
		}
		
	}

}