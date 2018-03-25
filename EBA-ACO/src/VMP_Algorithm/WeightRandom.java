package VMP_Algorithm;


import java.util.HashMap;  
import java.util.Map;  

class WeightRandom {  
    java.util.Random r = new java.util.Random();  
  
    public int getWeightRandom(double [] weightArrays,double weightSum) {  
        
        double stepWeightSum = 0;  
        for (int i = 0; i < weightArrays.length; i++) {  
            stepWeightSum += weightArrays[i];  
            if (Math.random() <= stepWeightSum/weightSum) {  
                //System.out.println(i);  
                return i;  
            }  
        }  
        
        for(double weight: weightArrays)
        {
        	System.out.print(weight + ",");
        }
        
        System.out.println(); 
        System.out.println("weightSum£º" + weightSum); 
        System.out.println("³ö´íÎóÁË");  
        return -1;  
    }     
}  