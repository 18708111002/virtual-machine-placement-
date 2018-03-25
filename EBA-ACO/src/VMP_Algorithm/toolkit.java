package VMP_Algorithm;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;


class GetFoldFileNames {

    public static List<String> getFileName(String path) {
       
        File f = new File(path);
        if (!f.exists()) {
            System.out.println(path + " not exists");
            return null;
        }
        List<String> fileList = new ArrayList<String>();

        File fa[] = f.listFiles();
        for (int i = 0; i < fa.length; i++) {
            File fs = fa[i];
            if (fs.isDirectory()) {
                System.out.println(fs.getName() + " [目录]");
            } else {
//                System.out.println(fs.getName());
                fileList.add(fs.getName());
            }
        }
        
        return fileList;
    }
}

class RequirementWR {

	public static void generateRequirement(double[][] arr, int n, int m,String fileName) throws IOException
	{
//		File file = new File("f:\\array.txt");  //存放数组数据的文件
		
		File file = new File(fileName);  //存放数组数据的文件
		 
		FileWriter out = new FileWriter(file);  //文件写入流
		 
		  //将数组中的数据写入到文件中。每行各数据之间TAB间隔
		for(int i = 0; i < n; i++){
			for(int j = 0; j < m;j++){
  
			   out.write(arr[i][j] + "\t");
		   }
		   out.write("\r\n");
		  }
		  out.close();
	}
	
	public static void readRequirement(double[][] arr,String fileName) throws IOException
	{
		BufferedReader in = new BufferedReader(new FileReader(fileName));  //
		String line;  //一行数据
		int row=0;
		
		  //逐行读取，并将每个数组放入到数组中
		 while((line = in.readLine()) != null){
			 
			 String[] temp = line.split("\t"); 
			 
			 for(int j = 0; j < temp.length; j++){
				 
				 arr[row][j] = Double.parseDouble(temp[j]);
		   }
			 row++;
		  }
		 in.close();
	}
	
	public static void generateRequirement(int[] arr, int n,int lineNum,String fileName) throws IOException
	{
	
		File file = new File(fileName);  
		FileWriter out = new FileWriter(file);  
		 
		  //将数组中的数据写入到文件中。每行各数据之间TAB间隔
		for(int i = 0; i < n; i++){
			
			out.write(arr[i] + "\t");
		   
			if((i + 1) % lineNum == 0)
				out.write("\r\n");
		  }
		  out.close();
	}
	
	public static void readRequirement(int[] arr,String fileName) throws IOException
	{
		BufferedReader in = new BufferedReader(new FileReader(fileName));  //
		String line;  //一行数据
		int count = 0;
		
		  //逐行读取，并将每个数组放入到数组中
		 while((line = in.readLine()) != null){
			 
			 String[] temp = line.split("\t"); 
			 
			 for(int j = 0; j < temp.length; j++){
				 
				 arr[count++] = Integer.parseInt(temp[j]);
		   }
			 
		  }
		 in.close();
	}
}

class RankIndex {

    public static List<Double> minN(List<Double> lst,int n,int []index) {
        if (lst.size() <= n)
            return lst;
        Double a = lst.remove(lst.size() - 1); 
        int aIndex = lst.size();
        int tIndex;
        
        List<Double> b = minN(lst,n,index);
        //System.out.println(b);
        for (int i = 0; i < b.size(); i++) {
            Double t = b.get(i);
            if (a < t) {
                //System.out.println(a + " : " + t);
                lst.set(i, a); 
                tIndex = index[i];
                index[i] = aIndex;
                aIndex = tIndex;   
                a = t;
            }
        }
        return b;
    }
}

class SortByCPU implements Comparator 
{

	@Override
	public int compare(Object arg0, Object arg1) {
		// TODO Auto-generated method stub
		vmRequire vmR1 = (vmRequire)arg0;
		vmRequire vmR2 = (vmRequire)arg1;
		
		return vmR2.cpu -vmR1.cpu;		
	}
	
}

class SortByCPU2 implements Comparator 
{

	@Override
	public int compare(Object arg0, Object arg1) {
		// TODO Auto-generated method stub
		vmRequire vmR1 = (vmRequire)arg0;
		vmRequire vmR2 = (vmRequire)arg1;
		
		return vmR1.cpu -vmR2.cpu;		
	}
	
}

class SortByMem implements Comparator 
{

	@Override
	public int compare(Object arg0, Object arg1) {
		// TODO Auto-generated method stub
		vmRequire vmR1 = (vmRequire)arg0;
		vmRequire vmR2 = (vmRequire)arg1;		
		return vmR2.mem -vmR1.mem;		
	}
	
}

class SortByTraffic implements Comparator 
{
	@Override
	public int compare(Object arg0, Object arg1) {
		// TODO Auto-generated method stub
		VMtraffic vmT1 = (VMtraffic)arg0;
		VMtraffic vmT2 = (VMtraffic)arg1;		
		
		 if((vmT1.traffic - vmT2.traffic) > 0)   
             return -1;  
         else if((vmT1.traffic - vmT2.traffic) < 0)  
             return 1;  
         else return 0; 				
	}	
}

class VMtraffic
{
	public double traffic;
	public int vm1;
	public int vm2;
	
	VMtraffic(int vm1,int vm2,double traffic){
		this.vm1 = vm1;
		this.vm2 = vm2;
		this.traffic = traffic;
	}
}
class vmRequire
{
	public int cpu;
	public int mem;
	public int id;
	
	vmRequire(int cpu,int mem,int id)
	{
		this.cpu = cpu;
		this.mem = mem;
		this.id  = id;
	}
}


class TestResultWR {

	public static void generateResult(String name,double[] arr, int n,int lineNum,String fileName) throws IOException
	{
		
		java.text.DecimalFormat df = new java.text.DecimalFormat("0.00");  
		File file = new File(fileName);  
		FileWriter out = new FileWriter(file,true);  
		out.write("\r\n");
		out.write(name + ":\r\n");
		
		int count = 0;
		 
		  //将数组中的数据写入到文件中。每行各数据之间TAB间隔
		for(int i = 0; i < n; i++){
			
			out.write(++count + " " + name + ": " +  df.format(arr[i]) + "\t");
		   
			if((i + 1) % lineNum == 0)
				out.write("\r\n");
		  }
		  out.close();
	}
	
	public static double getAverage(double [] arr,int num){
		
        double sum = 0;
        
        for(int i = 0;i < num;i++){
            sum += arr[i];
        }
        return (double)(sum / num);
    }

	public static double getStandardDevition(double [] array,int num,double average){
        double sum = 0;
        for(int i = 0;i < num;i++){
            sum +=((double)array[i] -average) * ((double)array[i] -average);
        }
        return Math.sqrt((sum / (num - 1)));
    }
	
}