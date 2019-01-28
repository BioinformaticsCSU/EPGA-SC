//package program;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

class CommonClass {
	//Get File Lines.
	public static int getFileLines(String ReadSetPath) throws IOException {
		int line = 0;
		String encoding = "utf-8";
		File file = new File(ReadSetPath);
		if (file.isFile() && file.exists()) {
			InputStreamReader read = new InputStreamReader(new FileInputStream(file), encoding);
			BufferedReader bufferedReader = new BufferedReader(read);
			while ((bufferedReader.readLine()) != null) {
				line++;
			}
			bufferedReader.close();
		}
		return line;
	}
	//File file to array.
	public static int FileToArray(String FileSetPath, String[] FileSetArray, String FilterChar) {
		int count = 0;
		String readtemp="";
		try {
			String encoding = "utf-8";
			File file = new File(FileSetPath);
			if (file.isFile() && file.exists()) {
				InputStreamReader read = new InputStreamReader(new FileInputStream(file), encoding);
				BufferedReader bufferedReader = new BufferedReader(read);
				while ((readtemp=bufferedReader.readLine())!= null) {
					if(!readtemp.startsWith(FilterChar) && readtemp.length()>5){
						FileSetArray[count++]=readtemp;
					}
				}
				bufferedReader.close();
			} else {
				System.out.println("File is not exist!");
			}
		} catch (Exception e) {
			System.out.println("Error liaoxingyu");
			e.printStackTrace();
		}
		return count;
	}
	//SAN File file to array.
	public static int SamFileToArray(String FileSetPath, String[] FileSetArray) {
		int count = 0;
		String readtemp="";
		try {
			String encoding = "utf-8";
			File file = new File(FileSetPath);
			if (file.isFile() && file.exists()) {
				InputStreamReader read = new InputStreamReader(new FileInputStream(file), encoding);
				BufferedReader bufferedReader = new BufferedReader(read);
				while ((readtemp=bufferedReader.readLine())!= null) {
					if(readtemp.charAt(0)!='@'){
				    	String [] Splitp = readtemp.split("\t|\\s+");
				    	if(!Splitp[1].equals("4")){
						   FileSetArray[count++]=readtemp;
				    	}
					}
				}
				bufferedReader.close();
			} else {
				System.out.println("File is not exist!");
			}
		} catch (Exception e) {
			System.out.println("Error liaoxingyu");
			e.printStackTrace();
		}
		return count;
	}
	//Estimate k-mer depth.
	public static double EstimatingKmerDepth(String gceFile,String WritePath) throws IOException {
	    double PeakValue=0;
		String encoding = "utf-8";
		String readtemp="";
		File file = new File(gceFile);
		if (file.isFile() && file.exists()) {
			InputStreamReader read = new InputStreamReader(new FileInputStream(file), encoding);
			BufferedReader bufferedReader = new BufferedReader(read);
			while ((readtemp=bufferedReader.readLine()) != null) {
				if(readtemp.equals("Final estimation table:")){
					bufferedReader.readLine();
					String GetString=bufferedReader.readLine();
					String [] SplitLine = GetString.split("\t|\\s+");
					PeakValue= Double.valueOf(SplitLine[4]);
					break;
				}
			}
			bufferedReader.close();
		}
		//Write.
		FileWriter writer_fasta= new FileWriter(WritePath+"/Estimated_kmerdepth.txt",true);
		writer_fasta.write("The peak value of k-mer histogram="+PeakValue+"\n");
		writer_fasta.close();
		return PeakValue;
	}
	//k-mer classification.
	public static int KmerClassification(String KmerFilePath,String LowFreKmerArray[],double LowDepthThreshold) throws IOException {
		double PeakValue=0;
		int CountLowFreKmer=0;
		String encoding = "utf-8";
		String readtemp="";
		File file = new File(KmerFilePath);
		if (file.isFile() && file.exists()) {
			InputStreamReader read = new InputStreamReader(new FileInputStream(file), encoding);
			BufferedReader bufferedReader = new BufferedReader(read);
			while ((readtemp=bufferedReader.readLine()) != null) {
				String [] SplitLine = readtemp.split("\t|\\s+");
				PeakValue= Double.valueOf(SplitLine[1]);
				if(PeakValue<LowDepthThreshold){
					LowFreKmerArray[CountLowFreKmer++]=SplitLine[0];
				}
			}
			bufferedReader.close();
		}
		return CountLowFreKmer;
	}
    //Generate FASTA From FASTQ.
    public static void GenerateFastaFromFastqFiles(String ParentPath,String Fq_left,String Fq_right,String DataName) {
		//TODO Auto-generated method stub
		try {
			  Process p_fq2fa_left=null;
			  Process p_fq2fa_right=null;
			  Runtime r_fq2fa_left=Runtime.getRuntime();
			  Runtime r_fq2fa_right=Runtime.getRuntime();
		      try{
			  	    String[] cmd_fq2fa_left = { "sh", "-c", ParentPath+"/tools/fq2fa "+Fq_left+" "+ParentPath+"/ReadFiles/"+DataName+"_left.fa"};
			  	    String[] cmd_fq2fa_right = { "sh", "-c", ParentPath+"/tools/fq2fa "+Fq_right+" "+ParentPath+"/ReadFiles/"+DataName+"_right.fa"};
			  	    p_fq2fa_left=r_fq2fa_left.exec(cmd_fq2fa_left);
			  	    p_fq2fa_left.waitFor();
			  	    p_fq2fa_right=r_fq2fa_right.exec(cmd_fq2fa_right);
			  	    p_fq2fa_right.waitFor();
		      }
		      catch(Exception e){
		    	    System.out.println("Step01 Error:"+e.getMessage());
		    	    e.printStackTrace();
		      }
		      //Merge FASTA.      
			  Process p_Merge=null;
			  Runtime r_Merge=Runtime.getRuntime();
		      try{
			  	    String[] cmd_gce = { "sh", "-c", ParentPath+"/tools/shuffleSequences_fasta.pl "+ParentPath+"/ReadFiles/"+DataName+"_left.fa"+" "+ParentPath+"/ReadFiles/"+DataName+"_right.fa"+" "+ParentPath+"/ReadFiles/"+DataName+".fa"};
			  	    p_Merge=r_Merge.exec(cmd_gce);
					p_Merge.waitFor();
		      }
		      catch(Exception e){
		    	    System.out.println("Step01 Error:"+e.getMessage());
		    	    e.printStackTrace();
		      }
		      
		} catch (Exception e) {
			System.out.println("Error liaoxingyu");
			e.printStackTrace();
		}
    }
    //Delete files.
    public static void deleteFile(File file){
        if (file.isFile()){
            file.delete();
        }else{
            String[] childFilePath = file.list();
            for (String path:childFilePath){
                File childFile= new File(file.getAbsoluteFile()+"/"+path);
                deleteFile(childFile);
            }
            System.out.println(file.getAbsoluteFile());
            file.delete();
        }
    }
    //Clean files.
    public static void delAllFile(File path) {
        if (!path.exists() || !path.isDirectory()) {
            return ;
        }
        String[] tmpList = path.list();
        if (tmpList != null) {
            for (String aTempList : tmpList) {
                File tmpFile = new File(path, aTempList);
                if (tmpFile.isFile() && (tmpFile.getName().endsWith(".gce")||tmpFile.getName().endsWith(".error")||tmpFile.getName().endsWith(".fa")||tmpFile.getName().endsWith(".h5")||tmpFile.getName().endsWith(".SAMFile")||tmpFile.getName().endsWith(".gz")||tmpFile.getName().endsWith(".depth")||tmpFile.getName().endsWith(".fastg")||tmpFile.getName().endsWith(".yaml")||tmpFile.getName().endsWith(".info")||tmpFile.getName().endsWith(".gfa")||tmpFile.getName().endsWith(".paths")||tmpFile.getName().endsWith(".fasta")||tmpFile.getName().endsWith(".bt2")||tmpFile.getName().endsWith(".fa")||tmpFile.getName().endsWith(".sam")||tmpFile.getName().endsWith(".bam"))||tmpFile.getName().endsWith(".bai")) {
                    tmpFile.delete();
                } else if (tmpFile.isDirectory()) {
                    delAllFile(tmpFile);
                }
            }
        }
    }
    //Clean special files.
    public static void delSpecialFile(File path, String Prefix, String suffix) {
        if (!path.exists() || !path.isDirectory()) {
            return ;
        }
        String[] tmpList = path.list();
        if (tmpList != null) {
            for (String aTempList : tmpList) {
                File tmpFile = new File(path, aTempList);
                if (tmpFile.isFile() && (tmpFile.getName().endsWith(suffix)||tmpFile.getName().startsWith(Prefix))) {
                    tmpFile.delete();
                }
            }
        }
    } 
	//Change lines of CONTIGS file.
	public static void MergeContigMultiLines(String ContigPath,String WritePath) throws IOException {
		int LineCount = 0;
		int LineNum = 0;
		String encoding = "utf-8";
		String readtemp = "";
		String ReadTemp = "";
		try {
			File file = new File(ContigPath);
			if (file.isFile() && file.exists()) {
				InputStreamReader read = new InputStreamReader(new FileInputStream(file), encoding); 
				BufferedReader bufferedReader = new BufferedReader(read);
				while ((readtemp = bufferedReader.readLine()) != null) {
					if (readtemp.charAt(0) == '>' && LineCount == 0) {
						ReadTemp=">NODE_" + (LineNum++) + "_Length" + "\n";
						LineCount++;
					} 
					else {
						if (readtemp.charAt(0) == '>'){
							FileWriter writer1= new FileWriter(WritePath,true);
							writer1.write(ReadTemp);
							writer1.close();
							ReadTemp="\n>NODE_"+(LineNum++)+"_Length"+"\n";
						} else {
							ReadTemp+= readtemp;
						}
					}
				}
				FileWriter writer1= new FileWriter(WritePath,true);
				writer1.write(ReadTemp);
				writer1.close();
				bufferedReader.close();
			} else {
				System.out.println("File is not exist!");
			}
		} catch (Exception e) {
			System.out.println("Error liaoxingyu");
			e.printStackTrace();
		}
	}
	//Reverse
	public static String reverse(String s) {
		int length = s.length();
		String reverse = "";
		for (int i = 0; i < length; i++) {
			if (s.charAt(i) == 'A') {
				reverse = "T" + reverse;
			} else if (s.charAt(i) == 'T') {
				reverse = "A" + reverse;
			} else if (s.charAt(i) == 'G') {
				reverse = "C" + reverse;
			} else if (s.charAt(i) == 'C') {
				reverse = "G" + reverse;
			} else {
				reverse = "N" + reverse;
			}
		}
		return reverse;
	}
	//Getting the longest common string. 
	public static String maxCommonstring(String strOne, String strTwo){
		if(strOne==null||strTwo==null){
			return null;
		}
		if(strOne.equals("")||strTwo.equals("")){
			return null;
		}
		String max = "";
		String min = "";
		if(strOne.length()<strTwo.length()){
			max = strTwo;
			min = strOne;
		} else{
			max = strTwo;
			min = strOne;
		}
		String current ="";
		for(int i=0;i<min.length();i++){
			for(int begin=0,end=min.length()-i;end<=min.length();begin++,end++){
				current=min.substring(begin, end);
				if(max.contains(current)){
					return current;
				}
			}
		}
		return null;
	}
	//Hash Function.
	public static int RSHash(String str) {
		int hash = 0;
		for (int i = 0; i < str.length(); i++) {
			hash = str.charAt(i) + (hash << 6) + (hash << 16) - hash;
		}
		return (hash & 0x7FFFFFFF);
	}
	//KMER File To HashTable.
	public static void KmerFileToHash(String KmerFilePath, String[] KmerHashTable, int SizeOfHash) throws IOException {
		String encoding = "utf-8";
		try {
			String readtemp = "";
			File file = new File(KmerFilePath);
			if (file.isFile() && file.exists()) {
				InputStreamReader read = new InputStreamReader(new FileInputStream(file), encoding);
				BufferedReader bufferedReader = new BufferedReader(read);
				while ((readtemp = bufferedReader.readLine()) != null) {
					if (readtemp.length()>10) {
						int hashCode = RSHash(readtemp) % SizeOfHash;
						if (KmerHashTable[hashCode]==null) {
							KmerHashTable[hashCode]=readtemp;
						} 
						else {
							int i = 1;
							while (KmerHashTable[(hashCode + i) % SizeOfHash]!= null) {
								i++;
							}
							if (KmerHashTable[(hashCode + i) % SizeOfHash]== null) {
								KmerHashTable[(hashCode + i) % SizeOfHash]=readtemp;
							}
					   }
				   }
				}
				bufferedReader.close();
			} else {
				System.out.println("File is not exist!");
			}
		} catch (Exception e) {
			System.out.println("Error liaoxingyu");
			e.printStackTrace();
		}
	}
	//Get element From HashTable.
	public static int getHashUnit(String KmerStr, String HashTable[], int SizeOfDSK) {
		int i = 1;
		int hashcode = RSHash(KmerStr) % SizeOfDSK;
		if (HashTable[hashcode]!= null) {
			if (HashTable[hashcode].equals(KmerStr)) {
				return hashcode;
			} else {
				while (HashTable[(hashcode + i) % SizeOfDSK]!= null) {
					if (HashTable[(hashcode + i) % SizeOfDSK].equals(KmerStr)) {
						return (hashcode + i) % SizeOfDSK;
					} else {
						i++;
					}
				}
			}
		}
		return -1;
	}
	//Set output format.
	public static String Changeformat(double value) {
		DecimalFormat df = new DecimalFormat("0.00");
		df.setRoundingMode(RoundingMode.HALF_UP);
		return df.format(value);
	}
	//Copy file.
	public static void copyFile(String oldPath, String newPath) throws IOException {
        File oldFile = new File(oldPath);
        File file = new File(newPath);
        FileInputStream in = new FileInputStream(oldFile);
        FileOutputStream out = new FileOutputStream(file);
        byte[] buffer=new byte[2097152];
        int readByte = 0;
        while((readByte = in.read(buffer)) != -1){
            out.write(buffer, 0, readByte);
        }
        in.close();
        out.close();
    }
	//FASTA file to array.
	public static int FastaToArray(String ReadSetPath, String[] ReadSetArray) {
		int count = 0;
		try {
			String encoding = "utf-8";
			File file = new File(ReadSetPath);
			if (file.isFile() && file.exists()) {
				InputStreamReader read = new InputStreamReader(new FileInputStream(file), encoding);
				BufferedReader bufferedReader = new BufferedReader(read);
				while ((bufferedReader.readLine()) != null) {
					ReadSetArray[count++] = bufferedReader.readLine();
				}
				bufferedReader.close();
			} else {
				System.out.println("File is not exist!");
			}
		} catch (Exception e) {
			System.out.println("Error liaoxingyu");
			e.printStackTrace();
		}
		return count;
	}
	//Generate File Array.
	public static int FileToArray(String FilePath, String[] FileArray) throws IOException {
		int ReadCount = 0;
		String encoding = "utf-8";
		try {
			String readtemp = "";
			File file = new File(FilePath);
			if (file.isFile() && file.exists()) {
				InputStreamReader read = new InputStreamReader(new FileInputStream(file), encoding);
				BufferedReader bufferedReader = new BufferedReader(read);
				while ((readtemp = bufferedReader.readLine()) != null) {
					FileArray[ReadCount++] = readtemp;
				}
				bufferedReader.close();
			} else {
				System.out.println("File is not exist!");
			}
		} catch (Exception e) {
			System.out.println("Error liaoxingyu");
			e.printStackTrace();
		}
		return ReadCount;
	}
	//break error in assemblies.
	public static void breakErrorPoints(int r,String home,String MUMmerFile,String EPGAcontigFile,String SPAdescontigFile,String FinalEPGAcontigPath, String FilterThreshold) throws IOException {
        int filterThreshold=Integer.parseInt(FilterThreshold);
		int SizeOfMUMmerFile=CommonClass.getFileLines(MUMmerFile);
        String MUMerArray[]=new String[SizeOfMUMmerFile];
        int RealSizeMUMmer=CommonClass.FileToArray(MUMmerFile,MUMerArray);
        System.out.print("[MUMmer_Num:"+RealSizeMUMmer);
        //Load EPGA.
		int SizeOfEPGAFile=CommonClass.getFileLines(EPGAcontigFile);
        String EPGAcontigArray[]=new String[SizeOfEPGAFile];
        int RealSizeEPGAcontig=CommonClass.FastaToArray(EPGAcontigFile,EPGAcontigArray);
        System.out.print(" EPGA _Num:"+RealSizeEPGAcontig);
        //Load Velvet.
		int SizeOfSPAdesFile=CommonClass.getFileLines(SPAdescontigFile);
        String SPAdescontigArray[]=new String[SizeOfSPAdesFile];
        int RealSizeSPAdescontig=CommonClass.FastaToArray(SPAdescontigFile,SPAdescontigArray);
        System.out.print(" SPAdes _Num:"+RealSizeSPAdescontig);
        //Process.
        int NumLines=0;
        String SaveTempArray[]=new String[RealSizeMUMmer];
        for(int w=4;w<RealSizeMUMmer;w++){
        	if(MUMerArray[w].charAt(0)!='#'){
	        	int CountSave=0;
	        	String MUMmerGroupRecords="";
	        	String [] SplitLine1 = MUMerArray[w].split("\t|\\s+");
	        	SaveTempArray[CountSave++]=MUMerArray[w];
	        	MUMmerGroupRecords=MUMerArray[w]+"\n";
	        	MUMerArray[w]="#"+MUMerArray[w];
	            for(int e=w+1;e<RealSizeMUMmer;e++){
	            	if(MUMerArray[e].charAt(0)!='#'){
		            	String [] SplitLine2 = MUMerArray[e].split("\t|\\s+");
		            	if(SplitLine1[11].equals(SplitLine2[11])){
		            		SaveTempArray[CountSave++]=MUMerArray[e];
		            		MUMmerGroupRecords+=MUMerArray[e]+"\n";
		            		MUMerArray[e]="#"+MUMerArray[e];
		            	}
	            	}
	            }
			    FileWriter writer1= new FileWriter(home+"/Alignments/MUMmerGroupRecords.fa",true);
                writer1.write(MUMmerGroupRecords+"-------------------------------------------\n");
                writer1.close();
                //Mask CONTIGS.
				String [] SplitLine3 = SaveTempArray[0].split("\t|\\s+");
	            String [] SplitLine4 = SplitLine3[11].split("--");
	            int EPGS_id=Integer.parseInt(SplitLine4[0]);
	            int SaveCoords[]=new int[EPGAcontigArray[EPGS_id].length()];
	            char SaveContigString[]=new char[EPGAcontigArray[EPGS_id].length()];
	            for(int h=0;h<EPGAcontigArray[EPGS_id].length();h++)
	            {
	            	SaveContigString[h]=EPGAcontigArray[EPGS_id].charAt(h);
	            }
                for(int f=0;f<CountSave;f++){
					String [] SplitLine5 = SaveTempArray[f].split("\t|\\s+");
		            int EPGA_start=Integer.parseInt(SplitLine5[0]);
		            int EPGA_end=Integer.parseInt(SplitLine5[1]);
		            int Cut_Start=EPGA_start-1;
		            int Cut_End=EPGA_end-1;
		            for(int h=Cut_Start;h<Cut_End;h++){
		            	SaveCoords[h]=1;
		            }
                }
                for(int c=0;c<EPGAcontigArray[EPGS_id].length();c++){
                	if(SaveCoords[c]==1){
                		SaveContigString[c]='N';
                	}
                }
	            String Str=new String(SaveContigString);
	            String [] SplitLine = Str.split("N");
	            for(int y=0;y<SplitLine.length;y++){
	            	if(SplitLine[y].length()>=filterThreshold)
	            	{
	   	               FileWriter writer2= new FileWriter(FinalEPGAcontigPath+"/contigs_corr.fa",true);
	                   writer2.write(">A"+(NumLines++)+"\n"+SplitLine[y]+"\n");
	                   writer2.close();
	            	}
	            }
        	}
        }
        for(int f=0;f<RealSizeSPAdescontig;f++)
        {
            FileWriter writer2= new FileWriter(FinalEPGAcontigPath+"/contigs_corr.fa",true);
            writer2.write(">M"+(NumLines++)+"\n"+SPAdescontigArray[f]+"\n");
            writer2.close();
        } 
	}
    //Add more common functions.
}
//Search Elements By Multiple Threads.
class MarkingReadByMultipleThreads implements Runnable{
	int SplitSize=0;
	String ReadArray[];
	String lowDepthKmerHashTable_A[];
	String lowDepthKmerHashTable_T[];
	String lowDepthKmerHashTable_G[];
	String lowDepthKmerHashTable_C[];
	int Start;
	int End;
	int ArraySize=0;
	int KmerSize=0;
	int Size_LowDepthHash=0;
	double threhold=0;
	String encoding="utf-8";
	//Construction.
	public MarkingReadByMultipleThreads(double Threhold,int kmerSize,int Count,int SizeOfLowDepthHash,String LowDepthKmerHashTable_A[],String LowDepthKmerHashTable_T[],String LowDepthKmerHashTable_G[],String LowDepthKmerHashTable_C[],String readArray[],int splitSize,int StartPosition,int Endposition)
	{
		KmerSize=kmerSize;
		Start=StartPosition;
		End=Endposition;
		SplitSize=splitSize;
		ReadArray=readArray;
		ArraySize=Count;
		threhold=Threhold;
		Size_LowDepthHash=SizeOfLowDepthHash;
		lowDepthKmerHashTable_A=LowDepthKmerHashTable_A;
		lowDepthKmerHashTable_T=LowDepthKmerHashTable_T;
		lowDepthKmerHashTable_G=LowDepthKmerHashTable_G;
		lowDepthKmerHashTable_C=LowDepthKmerHashTable_C;
	}
	@Override
	public void run() {
		for(int i=Start;i<=End;i++){
			if(ReadArray[i].charAt(0)!='#'){
				int CurrentLowFreKmerNum=0;
				int Num_Kmer=ReadArray[i].length()-KmerSize+1;
				String Rv_string=CommonClass.reverse(ReadArray[i]);
	    	    for(int c=0;c<=ReadArray[i].length()-KmerSize;c++){
	    	    	String FR_kmer=ReadArray[i].substring(c,c+KmerSize);
	    	    	String RF_kmer=Rv_string.substring(c,c+KmerSize);
	    	    	int Hash_Low_FR=0;
	    	    	int Hash_Low_RF=0;
	    	    	//FR.
	    	    	if(FR_kmer.charAt(0)=='A'){
	    	    		Hash_Low_FR=CommonClass.getHashUnit(FR_kmer,lowDepthKmerHashTable_A,Size_LowDepthHash);
	    	    	}
	    	    	if(FR_kmer.charAt(0)=='T'){
	    	    		Hash_Low_FR=CommonClass.getHashUnit(FR_kmer,lowDepthKmerHashTable_T,Size_LowDepthHash);
	    	    	}
	    	    	if(FR_kmer.charAt(0)=='G'){
	    	    		Hash_Low_FR=CommonClass.getHashUnit(FR_kmer,lowDepthKmerHashTable_G,Size_LowDepthHash);
	    	    	}
	    	    	if(FR_kmer.charAt(0)=='C'){
	    	    		Hash_Low_FR=CommonClass.getHashUnit(FR_kmer,lowDepthKmerHashTable_C,Size_LowDepthHash);
	    	    	}
	    	    	//RF.
	    	    	if(RF_kmer.charAt(0)=='A'){
	    	    		Hash_Low_RF=CommonClass.getHashUnit(RF_kmer,lowDepthKmerHashTable_A,Size_LowDepthHash);
	    	    	}
	    	    	if(RF_kmer.charAt(0)=='T'){
	    	    		Hash_Low_RF=CommonClass.getHashUnit(RF_kmer,lowDepthKmerHashTable_T,Size_LowDepthHash);
	    	    	}
	    	    	if(RF_kmer.charAt(0)=='G'){
	    	    		Hash_Low_RF=CommonClass.getHashUnit(RF_kmer,lowDepthKmerHashTable_G,Size_LowDepthHash);
	    	    	}
	    	    	if(RF_kmer.charAt(0)=='C'){
	    	    		Hash_Low_RF=CommonClass.getHashUnit(RF_kmer,lowDepthKmerHashTable_C,Size_LowDepthHash);
	    	    	}
	    	    	//Check.
	    	    	if(Hash_Low_FR!=-1||Hash_Low_RF!=-1){
	    	    		CurrentLowFreKmerNum++;
	    	    	}
	    	    	if(CurrentLowFreKmerNum>=(int)(Math.ceil((double)(threhold*Num_Kmer)))){
	    	    		ReadArray[i]="#"+ReadArray[i];
	    	    		break;
	    	    	}
	    	    } 
			}
	    }
	}
}
public class epga_sc {

	public static void main(String[] args) throws IOException, InterruptedException {
		// TODO Auto-generated method stub
		/***************************************************
		 ***************************************************
		 ***************************************************
		 * 'R':  The average length of read.
		 * 'c':  The k-mer size used in coverage estimate.
		 * 'k':  The k-mer size used in low depth reads assembly.
		 * 'K':  The k-mer size used in normal depth reads assembly.
		 * 't':  The number of threads.
		 * 'i':  The value of insert size.
		 * 's':  The standard devision of insert size.
		 * 'l':  The threshold of low depth k-mer (l*D=low depth threshold).
		 * 'r':  The threshold of the number of low depth k-mer in a low depth reads(r*(read_length-kmer_length+1)).
		 * 'q1': The first FASTQ file.
		 * 'q2': The second FASTQ file.
		 * 'o' : The path used to save final assemblies
		 ****************************************************
		 ****************************************************
		 ****************************************************/
		//definition of parameters.
		int    R = 100;
		int    c = 13;
		int    k = 21;
		int    K = 31;
		int    t = 64;
		int    i = 330;
		double s = 30;
		double l = 0.1;
		double r = 0.1;
		String q1 = "";
		String q2 = "";
		String o  = "";
		// get the value of each parameter.
		for (int x = 0; x < args.length; x += 2) {
			String headStr = args[x].substring(1, args[x].length());
			switch (headStr) {
			case "R": {
				R = Integer.parseInt(args[x + 1]);
			};break;
			case "c": {
				c = Integer.parseInt(args[x + 1]);
			};break;
			case "k": {
				k = Integer.parseInt(args[x + 1]);
			};break;
			case "K": {
				K = Integer.parseInt(args[x + 1]);
			};break;
			case "t": {
				t = Integer.parseInt(args[x + 1]);
			};break;
			case "i": {
				i = Integer.parseInt(args[x + 1]);
			};break;
			case "s": {
				s = Double.valueOf(args[x + 1]);
			};break;
			case "l": {
				l = Double.valueOf(args[x + 1]);
			};break;
			case "r": {
				r = Double.valueOf(args[x + 1]);
			};break;
			case "q1": {
				q1 = args[x + 1];
			};break;
			case "q2": {
				q2 = args[x + 1];
			};break;
			case "o": {
				o = args[x + 1];
			};break;
		 }
	  }
	  if(!q1.equals("")&&!q2.equals("")&&!o.equals(""))
	  {
          //Main Program begin.
	      System.out.println("-------------------------------------------------------------------------------------------");
	      System.out.println("Copyright (C) 2017 Jianxin Wang(jxwang@mail.csu.edu.cn), Xingyu Liao(liaoxingyu@csu.edu.cn)"+"\n"+"School of Information Science and Engineering"+"\n"+"Central South University"+"\n"+"ChangSha"+"CHINA, 410083"+"\n");		  
		  File directory = new File(".");
		  String ParentPath=directory.getCanonicalPath();
		  //Delete some required directories.
		  File CorrectedDirectoy = new File(ParentPath+"/Assembly");
		  if(CorrectedDirectoy.exists())
		  {
			  CommonClass.delAllFile(CorrectedDirectoy);
		  }
		  File AlignDirectoy = new File(ParentPath+"/Auxiliary_lib");
		  if(AlignDirectoy.exists())
		  {
			  CommonClass.delSpecialFile(AlignDirectoy,"DATA",".fq");
			  CommonClass.delSpecialFile(AlignDirectoy,"DATA",".aln");
		  }
		  File OutputDirectoy = new File(o);
		  if(OutputDirectoy.exists())
		  {
			  CommonClass.delSpecialFile(OutputDirectoy,"Contigs.",".fasta");
			  CommonClass.delSpecialFile(OutputDirectoy,"Scaffolds.",".fasta");
		  }
		  File QSDirectoy = new File(ParentPath+"/QSstatistics");
		  if(QSDirectoy.exists())
		  {
			  CommonClass.delAllFile(QSDirectoy);
		  }
		  File ReadFileDirectoy = new File(ParentPath+"/ReadFiles");
		  if(ReadFileDirectoy.exists())
		  {
			  CommonClass.delAllFile(ReadFileDirectoy);
		  }
		  File KmerFileDirectoy = new File(ParentPath+"/KmerFiles");
		  if(KmerFileDirectoy.exists())
		  {
			  CommonClass.delAllFile(KmerFileDirectoy);
		  }
		  File AlignFileDirectoy = new File(ParentPath+"/Alignments");
		  if(AlignFileDirectoy.exists())
		  {
			  CommonClass.delAllFile(AlignFileDirectoy);
		  }
		  File SSPACEFileDirectoy = new File(ParentPath+"/scaffolds_sspace");
		  if(SSPACEFileDirectoy.exists())
		  {
			  CommonClass.delAllFile(SSPACEFileDirectoy);
		  }
		  File LogFileDirectoy = new File(ParentPath+"/Log");
		  if(LogFileDirectoy.exists())
		  {
			  CommonClass.delAllFile(LogFileDirectoy);
		  }
		  //Create the required directories.
		  //k-mer file path.
	      String kmerpath = ParentPath +"/KmerFiles";
	      File kmerFileDirectory =new File(kmerpath);  
	      if  (!kmerFileDirectory .isDirectory()){ 
	    	  kmerFileDirectory.mkdir();
	      }
		  //Read file path.
	      String ReadFilepath = ParentPath +"/ReadFiles";
	      File ReadFileDirectory =new File(ReadFilepath);  
	      if  (!ReadFileDirectory.isDirectory()){ 
	    	  ReadFileDirectory.mkdir();
	      }
		  //QS file path.
	      String QSFilepath = ParentPath +"/QSstatistics";
	      File QSFileDirectory =new File(QSFilepath);  
	      if  (!QSFileDirectory.isDirectory()){ 
	    	  QSFileDirectory.mkdir();
	      }
	      //Log file path.
	      String Logpath = ParentPath +"/Log";
	      File LogFileDirectory =new File(Logpath);  
	      if  (!LogFileDirectory.isDirectory()){ 
	    	  LogFileDirectory.mkdir();
	      }
	      //Alignment file path.
	      String Alignfilepath = ParentPath +"/Alignments";
	      File AlignFileDirectory =new File(Alignfilepath);  
	      if  (!AlignFileDirectory.isDirectory()){ 
	    	  AlignFileDirectory.mkdir();
	      }
	      //SSPACE file path.
	      String SSPACEfilepath = ParentPath +"/scaffolds_sspace";
	      File SSPACEFileDirectory =new File(SSPACEfilepath);  
	      if  (!SSPACEFileDirectory.isDirectory()){ 
	    	  SSPACEFileDirectory.mkdir();
	      }
	      //Left FASTQ file.
	      String LeftFastapath = ParentPath +"/ReadFiles/read_left.fa";
	      File LeftFastaFile =new File(LeftFastapath);  
	      if  (!LeftFastaFile.isFile()){ 
	    	  LeftFastaFile.createNewFile();
	      }
	      //Right FASTQ file.
	      String RightFastapath = ParentPath +"/ReadFiles/read_right.fa";
	      File RightFastaFile =new File(RightFastapath);  
	      if  (!RightFastaFile.isFile()){ 
	    	  RightFastaFile.createNewFile();
	      }
	      //FASTA file path.
	      String Fastapath = ParentPath +"/ReadFiles/read.fa";
	      File FastaFile =new File(Fastapath);  
	      if  (!FastaFile.isFile()){ 
	    	  FastaFile.createNewFile();
	      }
	      //Assembly save path.
	      String Assemblypath = ParentPath +"/Assembly";
	      File AssemblyDirectory =new File(Assemblypath);  
	      if  (!AssemblyDirectory.isDirectory()){ 
	    	  AssemblyDirectory.mkdir();
	      }
	      //Assembly save path.
	      String SPAdesAssemblypath = ParentPath +"/Assembly/SPAdes";
	      File SPAdesAssemblyDirectory =new File(SPAdesAssemblypath);  
	      if  (!SPAdesAssemblyDirectory.isDirectory()){ 
	    	  SPAdesAssemblyDirectory.mkdir();
	      }
	      //Assembly save path.
	      String EPGA_SC_Assemblypath = ParentPath +"/Assembly/EPGA-SC";
	      File EPGA_SC_AssemblyDirectory =new File(EPGA_SC_Assemblypath);  
	      if  (!EPGA_SC_AssemblyDirectory.isDirectory()){ 
	    	  EPGA_SC_AssemblyDirectory.mkdir();
	      }
	      //Low depth reads assembly save path.
	      String EPGA_SC_Low_Assemblypath = ParentPath +"/Assembly/EPGA-SC/LowDepthReads";
	      File EPGA_SC_Low_AssemblyDirectory =new File(EPGA_SC_Low_Assemblypath);  
	      if  (!EPGA_SC_Low_AssemblyDirectory.isDirectory()){ 
	    	  EPGA_SC_Low_AssemblyDirectory.mkdir();
	      }
	      //Normal depth reads assembly save path.
	      String EPGA_SC_Normal_Assemblypath = ParentPath +"/Assembly/EPGA-SC/NormalDepthReads";
	      File EPGA_SC_Normal_AssemblyDirectory =new File(EPGA_SC_Normal_Assemblypath);  
	      if  (!EPGA_SC_Normal_AssemblyDirectory.isDirectory()){ 
	    	  EPGA_SC_Normal_AssemblyDirectory.mkdir();
	      }
	      //Final assembly save path.
	      String EPGA_SC_Final_Assemblypath = ParentPath +"/Assembly/EPGA-SC/FinalAssembly";
	      File EPGA_SC_Final_AssemblyDirectory =new File(EPGA_SC_Final_Assemblypath);  
	      if  (!EPGA_SC_Final_AssemblyDirectory.isDirectory()){ 
	    	  EPGA_SC_Final_AssemblyDirectory.mkdir();
	      }
	      //Auxiliary library save path.
	      String Auxiliarypath = ParentPath +"/Auxiliary_lib";
	      File AuxiliaryDirectory =new File(Auxiliarypath);  
	      if  (!AuxiliaryDirectory.isDirectory()){ 
	    	  AuxiliaryDirectory.mkdir();
	      }
	      //QS file save path.
	      String QSstatisticspath = ParentPath +"/QSstatistics";
	      File QSstatisticsDirectory =new File(QSstatisticspath);  
	      if  (!QSstatisticsDirectory.isDirectory()){ 
	    	  QSstatisticsDirectory.mkdir();
	      }
		  //Save Final Assemblies.
	      String FinalAssemblyPath = o;
	      File FinalAssemblyDirectory =new File(FinalAssemblyPath);  
	      if  (!FinalAssemblyDirectory.isDirectory()){ 
	    	  FinalAssemblyDirectory.mkdir();
	      }
	      //Step00-PRE-treatment.
	      System.out.print("Step01: Parameters configuration");
	      System.out.println(" [ R="+R+" c="+c+" k="+k+" K="+K+" t="+t+" i="+i+" s="+s+" l="+l+" r="+r+" ]");
		  String [] q1_split = q1.split("/");
		  String DataName1="";
		  if(q1_split[q1_split.length-1].endsWith(".fastq")){
			  String [] q1_name = q1_split[q1_split.length-1].split(".fastq");
		      DataName1=q1_name[0];
		  }
	      if(q1_split[q1_split.length-1].endsWith(".fq")){
		      String [] q1_name = q1_split[q1_split.length-1].split(".fq");
		      DataName1=q1_name[0];
		  }
		  String [] q2_split = q2.split("/");
		  String DataName2="";
		  if(q2_split[q2_split.length-1].endsWith(".fastq")){
			  String [] q2_name = q2_split[q2_split.length-1].split(".fastq");
		      DataName2=q2_name[0];
		  }
	      if(q2_split[q2_split.length-1].endsWith(".fq")){
		      String [] q2_name = q2_split[q2_split.length-1].split(".fq");
		      DataName2=q2_name[0];
		  }
		  File ReadFiles=new File(ParentPath+"/readfile/read.fa");
	      if(ReadFiles.exists())
	      {
	    	  CommonClass.deleteFile(ReadFiles);
	      }
	      File ReadMergeFiles=new File(ParentPath+"/ReadFiles/read.fa");
	      if(ReadMergeFiles.exists())
	      {
	    	  CommonClass.deleteFile(ReadMergeFiles);
	      }
	      File ReadFiles_left=new File(ParentPath+"/ReadFiles/read_left.fa");
	      if(ReadFiles_left.exists())
	      {
	    	  CommonClass.deleteFile(ReadFiles_left);
	      }
	      File ReadFiles_right=new File(ParentPath+"/ReadFiles/read_right.fa");
	      if(ReadFiles_right.exists())
	      {
	          CommonClass.deleteFile(ReadFiles_right);
	      }
	      File KmerFiles=new File(ParentPath+"/read.txt");
	      if(KmerFiles.exists())
	      {
	    	  CommonClass.deleteFile(KmerFiles);
	      }
	      File Kmerh5Files=new File(ParentPath+"/read.h5");
	      if(Kmerh5Files.exists())
	      {
	    	  CommonClass.deleteFile(Kmerh5Files);
	      }
	      File KmerHisFiles=new File(ParentPath+"/KmerFiles/kmer_histogram.txt");
	      if(KmerHisFiles.exists())
	      {
	    	  CommonClass.deleteFile(KmerHisFiles);
	      }
	      File HighFreKmerFiles=new File(ParentPath+"/KmerFiles/HighFrequencyKmer.fa");
	      if(HighFreKmerFiles.exists())
	      {
	    	  CommonClass.deleteFile(HighFreKmerFiles);
	      }
	      File GCEFiles_Error=new File(ParentPath+"/KmerFiles/gce.error");
	      if(GCEFiles_Error.exists())
	      {
	    	  CommonClass.deleteFile(GCEFiles_Error);
	      }
	      File GCEFiles_Table=new File(ParentPath+"/KmerFiles/gce.table");
	      if(GCEFiles_Table.exists())
	      {
	    	  CommonClass.deleteFile(GCEFiles_Table);
	      }
		  File ESTFiles=new File(ParentPath+"/Estimated_kmerdepth.txt");
	      if(ESTFiles.exists())
	      {
	    	  CommonClass.deleteFile(ESTFiles);
	      }
		  File SpadesScaffold=new File(ParentPath+"/Assembly/SPAdes/scaffolds.fasta");
	      if(SpadesScaffold.exists())
	      {
	    	  CommonClass.deleteFile(SpadesScaffold);
	      }
		  File ScaffoldChange=new File(ParentPath+"/Assembly/SPAdes/ChangeLine_Scaffolds.fa");
	      if(ScaffoldChange.exists())
	      {
	    	  CommonClass.deleteFile(ScaffoldChange);
	      }
	      //Step02-start.
  	      File read_txt=new File(ParentPath);
  	      CommonClass.delSpecialFile(read_txt,"read",".h5");
  	      CommonClass.delSpecialFile(read_txt,"read",".txt");
  	      CommonClass.delSpecialFile(read_txt,"Estimated_kmerdepth",".txt");
	      System.out.print("Step02: Data format conversion");
	      long startTime1 = System.currentTimeMillis(); 
	      Runtime r1 = Runtime.getRuntime();
	      long startMem1 = r1.freeMemory();
	      CommonClass.GenerateFastaFromFastqFiles(ParentPath,q1,q2,"read");
	      long orz1 = Math.abs(startMem1 - r1.freeMemory());
	      long endTime1 = System.currentTimeMillis();
	      System.out.println(" [Time consumption:"+(endTime1-startTime1)+"ms. Memory consumption:"+(double)orz1/1000000000+"GB] ");
	      //Step03-start.
	      System.out.print("Step03: K-mer frequency statistics");
	      long startTime2 = System.currentTimeMillis(); 
	      Runtime r2 = Runtime.getRuntime();
	      long startMem2 = r2.freeMemory();
	      //Run DSK.
		  Process p_dsk1=null;
	      Runtime r_dsk1=Runtime.getRuntime();
	      try{
	    	    String cmd1=ParentPath+"/tools/dsk -file "+ParentPath+"/ReadFiles/read.fa -kmer-size "+c+" -abundance-min 1";
	    	    p_dsk1=r_dsk1.exec(cmd1);
				p_dsk1.waitFor();
	      }
	      catch(Exception e){
	    	    System.out.println("Step03 Error:"+e.getMessage());
	    	    e.printStackTrace();
	      }
		  Process p_dsk2=null;
		  Runtime r_dsk2=Runtime.getRuntime();
	      try{
		  	    String cmd2=ParentPath+"/tools/dsk2ascii -file read.h5 -out read.txt";
	    	    p_dsk2=r_dsk2.exec(cmd2);
				p_dsk2.waitFor();
	      }
	      catch(Exception e){
	    	    System.out.println("Step03 Error:"+e.getMessage());
	    	    e.printStackTrace();
	      }
          //Counting
		  String encoding = "utf-8";
		  String ReadTemp="";
		  TreeMap<Integer, Integer> kmerFremap = new TreeMap<Integer, Integer>();
		  File file = new File(ParentPath+"/read.txt");
		  if (file.isFile() && file.exists()) {
				InputStreamReader read = new InputStreamReader(new FileInputStream(file), encoding);
				BufferedReader bufferedReader = new BufferedReader(read);
				while ((ReadTemp=bufferedReader.readLine())!= null) {
					String [] SplitLine = ReadTemp.split("\t|\\s+");
					int FreNum=Integer.parseInt(SplitLine[1]);
					kmerFremap.put(FreNum,kmerFremap.get(FreNum) == null ? 1 : kmerFremap.get(FreNum)+ 1);
				}
				bufferedReader.close();
		  }
		  //Write histogram.
	      Set<Integer> set1 = kmerFremap.keySet();
	      for (Object obj : set1) {
				FileWriter writer1= new FileWriter(ParentPath+"/KmerFiles/kmer_histogram.txt",true);
				writer1.write(obj+"\t"+kmerFremap.get(obj)+"\n");
				writer1.close();
	      }
	      //Free.
	      kmerFremap=null;
	      set1=null;
	      long orz2 = Math.abs(startMem2 - r2.freeMemory());
	      long endTime2 = System.currentTimeMillis();
	      System.out.println(" [Time consumption:"+(endTime2-startTime2)+"ms. Memory consumption:"+(double)orz2/1000000000+"GB] ");
	      //Step04-start.  
	      System.out.println("Step04: Estimating the average coverage of k-mers");
	      long startTime_tgce = System.currentTimeMillis(); 
	      Runtime r_tgce = Runtime.getRuntime();
	      long startMem_tgce = r_tgce.freeMemory();
		  System.out.println("=================Process Start==================");
		  Process p_gce=null;
		  Runtime r_gce=Runtime.getRuntime();
	      try{
		  	    String[] cmd_gce = { "sh", "-c", ParentPath+"/tools/gce -f "+ParentPath+"/KmerFiles/kmer_histogram.txt >"+ParentPath+"/KmerFiles/gce.table 2>"+ParentPath+"/KmerFiles/gce.error"};
		  	    p_gce=r_gce.exec(cmd_gce);
				p_gce.waitFor();
	      }
	      catch(Exception e){
	    	    System.out.println("Step04 Error:"+e.getMessage());
	    	    e.printStackTrace();
	      }
	      double KmerDepth=CommonClass.EstimatingKmerDepth(ParentPath+"/KmerFiles/gce.error",ParentPath);
		  double ReadDepth=(KmerDepth*R)/(R-k+1);
		  System.out.print("[K-mer coverage:"+KmerDepth+" Read Coverage:"+ReadDepth);
		  long orz_tgce = Math.abs(startMem_tgce - r_tgce.freeMemory());
	      long endTime_tgce = System.currentTimeMillis();
	      System.out.println(" Time consumption:"+(endTime_tgce-startTime_tgce)+"ms. Memory consumption:"+(double)orz_tgce/1000000000+"GB] ");
	      System.out.println("=================Process  End==================");
	      //Step05-Assembly.
	      System.out.print("Step05: SPAdes assembly");
		  long startTime_assembly = System.currentTimeMillis(); 
	      Runtime r_assembly = Runtime.getRuntime();
	      long startMem_assembly = r_assembly.freeMemory();
	      //Assembly start.
		  Process p_assemblyRun=null;
		  Runtime r_assemblyRun=Runtime.getRuntime();
	      try{
	    	    String[] cmd_assemblyRun = { "sh", "-c", ParentPath+"/tools/SPAdes-3.10.1-Linux/bin/spades.py -t 64 --pe1-1 "+q1+" --pe1-2 "+q2+" --sc --careful -o "+ParentPath+"/Assembly/SPAdes/ > "+ParentPath+"/Log/SPAdeslog.out"};
	    	    p_assemblyRun=r_assemblyRun.exec(cmd_assemblyRun);
				p_assemblyRun.waitFor();
	      }
	      catch(Exception e){
	    	    System.out.println("Step05 Error:"+e.getMessage());
	    	    e.printStackTrace();
	      }
		  //Assembly end.
		  long orz_assembly = Math.abs(startMem_assembly - r_assembly.freeMemory());
	      long endTime_assembly = System.currentTimeMillis();
	      System.out.println(" [Time consumption:"+(endTime_assembly-startTime_assembly)+"ms. Memory consumption:"+(double)orz_assembly/1000000000+"GB] ");
	      //Step06-Change lines of the Scaffolds.
	      System.out.print("Step06: Change lines of the Scaffolds");
	      long startTime_scaff = System.currentTimeMillis(); 
	      Runtime r_scaff = Runtime.getRuntime();
	      long startMem_scaff = r_scaff.freeMemory();
	      CommonClass.MergeContigMultiLines(ParentPath+"/Assembly/SPAdes/scaffolds.fasta",ParentPath+"/Assembly/SPAdes/scaffolds.ChangeLines.fasta");
	      long orz_scaff = Math.abs(startMem_scaff - r_scaff.freeMemory());
	      long endTime_scaff = System.currentTimeMillis();
	      System.out.println("[Time consumption:"+(endTime_scaff-startTime_scaff)+"ms. Memory consumption:"+(double)orz_scaff/1000000000+"GB] ");
	      //Step07-K-mer Classification.
	      System.out.print("Step07: k-mer classification");
	      long startTime_kc = System.currentTimeMillis(); 
	      Runtime r_kc = Runtime.getRuntime();
	      long startMem_kc = r_kc.freeMemory();
	      double LowFreThreshold=l*KmerDepth;
	      int LineOfDSKFile=CommonClass.getFileLines(ParentPath+"/read.txt");
	      String LowFreKmerArray[]=new String[LineOfDSKFile];
	      int SizeOfLowFreKmerFile=CommonClass.KmerClassification(ParentPath+"/read.txt",LowFreKmerArray,LowFreThreshold);
          for(int f=0;f<SizeOfLowFreKmerFile;f++){
        	    if(LowFreKmerArray[f].charAt(0)=='A')
        	    {
			       FileWriter writer1= new FileWriter(ParentPath+"/KmerFiles/LowFrequencyKmer_A.fa",true);
				   writer1.write(LowFreKmerArray[f]+"\n");
				   writer1.close();
        	    }
        	    if(LowFreKmerArray[f].charAt(0)=='T')
        	    {
			       FileWriter writer1= new FileWriter(ParentPath+"/KmerFiles/LowFrequencyKmer_T.fa",true);
				   writer1.write(LowFreKmerArray[f]+"\n");
				   writer1.close();
        	    }
        	    if(LowFreKmerArray[f].charAt(0)=='G')
        	    {
			       FileWriter writer1= new FileWriter(ParentPath+"/KmerFiles/LowFrequencyKmer_G.fa",true);
				   writer1.write(LowFreKmerArray[f]+"\n");
				   writer1.close();
        	    }
        	    if(LowFreKmerArray[f].charAt(0)=='C')
        	    {
			       FileWriter writer1= new FileWriter(ParentPath+"/KmerFiles/LowFrequencyKmer_C.fa",true);
				   writer1.write(LowFreKmerArray[f]+"\n");
				   writer1.close();
        	    }
			    FileWriter writer2= new FileWriter(ParentPath+"/KmerFiles/LowFrequencyKmer_All.fa",true);
				writer2.write(LowFreKmerArray[f]+"\n");
				writer2.close();
          }
	      long orz_kc = Math.abs(startMem_kc - r_kc.freeMemory());
	      long endTime_kc = System.currentTimeMillis();
	      System.out.println(" [Time consumption:"+(endTime_kc-startTime_kc)+"ms. Memory consumption:"+orz_kc/1000000000+"GB] ");
	      //Step08-Read Classification.
	      System.out.print("Step08: Read classification");
	      long startTime_rc=System.currentTimeMillis(); 
	      Runtime r_rc=Runtime.getRuntime();
	      long startMem_rc=r_kc.freeMemory();
	      //Left.
		  int ArraySize_left=CommonClass.getFileLines(ParentPath+"/ReadFiles/read_left.fa")/2;
		  String[] ReadSetArray_left=new String[ArraySize_left];
		  int scount_left=CommonClass.FileToArray(ParentPath+"/ReadFiles/read_left.fa",ReadSetArray_left,">");	
		  System.out.print("  [ArraySize_left:"+scount_left+" ");
		  //Right.
		  int ArraySize_right=CommonClass.getFileLines(ParentPath+"/ReadFiles/read_right.fa")/2;
		  String[] ReadSetArray_right=new String[ArraySize_right];
		  int scount_right=CommonClass.FileToArray(ParentPath+"/ReadFiles/read_right.fa",ReadSetArray_right,">");	
		  System.out.print(" ArraySize_right:"+scount_right+" ");
		  //LowFre HashTable
		  int SizeOfLowFreDSK=CommonClass.getFileLines(ParentPath +"/KmerFiles/LowFrequencyKmer_All.fa")+10000000;
		  System.out.print(" LowFreKmerSize:"+SizeOfLowFreDSK);
		  String LowFreHashTable_A[]=new String[SizeOfLowFreDSK];
		  String LowFreHashTable_T[]=new String[SizeOfLowFreDSK];
		  String LowFreHashTable_G[]=new String[SizeOfLowFreDSK];
		  String LowFreHashTable_C[]=new String[SizeOfLowFreDSK];
		  CommonClass.KmerFileToHash(ParentPath+"/KmerFiles/LowFrequencyKmer_A.fa",LowFreHashTable_A,SizeOfLowFreDSK);
		  CommonClass.KmerFileToHash(ParentPath+"/KmerFiles/LowFrequencyKmer_T.fa",LowFreHashTable_A,SizeOfLowFreDSK);
		  CommonClass.KmerFileToHash(ParentPath+"/KmerFiles/LowFrequencyKmer_G.fa",LowFreHashTable_A,SizeOfLowFreDSK);
		  CommonClass.KmerFileToHash(ParentPath+"/KmerFiles/LowFrequencyKmer_C.fa",LowFreHashTable_A,SizeOfLowFreDSK);
	      //Multiple Threads.
	      ExecutorService pool_1 = Executors.newFixedThreadPool(t);
	      int SplitSize = (int)(Math.ceil((double)scount_left/t));	    
		  for(int w=0;w<t;w++){
				if(w!=t-1){
					MarkingReadByMultipleThreads mt_l = new MarkingReadByMultipleThreads(r,c,scount_left,SizeOfLowFreDSK,LowFreHashTable_A,LowFreHashTable_T,LowFreHashTable_G,LowFreHashTable_C,ReadSetArray_left,SplitSize,w*SplitSize,(w+1)*SplitSize-1);
				    pool_1.execute(mt_l);
				}
				else{
					MarkingReadByMultipleThreads mt_l = new MarkingReadByMultipleThreads(r,c,scount_left,SizeOfLowFreDSK,LowFreHashTable_A,LowFreHashTable_T,LowFreHashTable_G,LowFreHashTable_C,ReadSetArray_left,SplitSize,w*SplitSize,scount_left-1);
					pool_1.execute(mt_l);
				}
		  }
		  //Waiting for all threads to complete.
		  pool_1.shutdown();
		  pool_1.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS); 
		  System.out.print(" Lib1-Over ");
		  //Multiple Threads.
	      ExecutorService pool_2 = Executors.newFixedThreadPool(t);
	      int SplitSize2 = (int)(Math.ceil((double)scount_right/t));	    
		  for(int w=0;w<t;w++){
				if(w!=t-1){
					MarkingReadByMultipleThreads mt_r = new MarkingReadByMultipleThreads(r,c,scount_right,SizeOfLowFreDSK,LowFreHashTable_A,LowFreHashTable_T,LowFreHashTable_G,LowFreHashTable_C,ReadSetArray_right,SplitSize2,w*SplitSize2,(w+1)*SplitSize2-1);
				    pool_2.execute(mt_r);
				}
				else{
					MarkingReadByMultipleThreads mt_r = new MarkingReadByMultipleThreads(r,c,scount_right,SizeOfLowFreDSK,LowFreHashTable_A,LowFreHashTable_T,LowFreHashTable_G,LowFreHashTable_C,ReadSetArray_right,SplitSize2,w*SplitSize2,scount_right-1);
					pool_2.execute(mt_r);
				}
		  }
		  //Waiting for all threads to complete.
		  pool_2.shutdown();
		  pool_2.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS); 
		  System.out.print(" Lib2-Over ");
          //OutPut low depth reads and normal depth reads.
		  int IndexLowDepthRead_left=0;
		  int IndexLowDepthRead_right=0;
		  int IndexNormalDepthRead_left=0;
		  int IndexNormalDepthRead_right=0;
		  String ReplaceString="N";
		  for(int g=1;g<R;g++){
			  ReplaceString+="N";
		  }
		  for(int z=0;z<scount_left;z++){
			   if(ReadSetArray_left[z].charAt(0)=='#' && ReadSetArray_right[z].charAt(0)=='#'){
				   FileWriter writer1= new FileWriter(ParentPath+"/ReadFiles/LowDepthReads_left.fa",true);
				   writer1.write(">Node_"+(IndexLowDepthRead_left++)+"\n"+ReadSetArray_left[z].substring(1,ReadSetArray_left[z].length())+"\n");
				   writer1.close();
				   FileWriter writer2= new FileWriter(ParentPath +"/ReadFiles/LowDepthReads_right.fa",true);
				   writer2.write(">Node_"+(IndexLowDepthRead_right++)+"\n"+ReadSetArray_right[z].substring(1,ReadSetArray_right[z].length())+"\n");
				   writer2.close();
			   }else{
				   if(ReadSetArray_left[z].charAt(0)!='#' && ReadSetArray_right[z].charAt(0)!='#'){
					   FileWriter writer1= new FileWriter(ParentPath+"/ReadFiles/NormalDepthReads_left.fa",true);
					   writer1.write(">Node_"+(IndexNormalDepthRead_left++)+"\n"+ReadSetArray_left[z].substring(0,ReadSetArray_left[z].length()-((int)(0.2*ReadSetArray_left[z].length())))+"\n");
					   writer1.close();
					   FileWriter writer2= new FileWriter(ParentPath +"/ReadFiles/NormalDepthReads_right.fa",true);
					   writer2.write(">Node_"+(IndexNormalDepthRead_right++)+"\n"+ReadSetArray_right[z].substring(0,ReadSetArray_right[z].length()-((int)(0.2*ReadSetArray_right[z].length())))+"\n");
					   writer2.close();
				   }
				   if(ReadSetArray_left[z].charAt(0)=='#' && ReadSetArray_right[z].charAt(0)!='#'){
					   FileWriter writer1= new FileWriter(ParentPath+"/ReadFiles/LowDepthReads_left.fa",true);
					   writer1.write(">Node_"+(IndexLowDepthRead_left++)+"\n"+ReadSetArray_left[z].substring(1,ReadSetArray_left[z].length())+"\n");
					   writer1.close();
					   FileWriter writer2= new FileWriter(ParentPath +"/ReadFiles/LowDepthReads_right.fa",true);
					   writer2.write(">Node_"+(IndexLowDepthRead_right++)+"\n"+ReadSetArray_right[z]+"\n");
					   writer2.close();
					   FileWriter writer3= new FileWriter(ParentPath+"/ReadFiles/NormalDepthReads_left.fa",true);
					   writer3.write(">Node_"+(IndexNormalDepthRead_left++)+"\n"+ReplaceString.substring(0,ReplaceString.length()-((int)(0.2*ReplaceString.length())))+"\n");
					   writer3.close();
					   FileWriter writer4= new FileWriter(ParentPath +"/ReadFiles/NormalDepthReads_right.fa",true);
					   writer4.write(">Node_"+(IndexNormalDepthRead_right++)+"\n"+ReadSetArray_right[z].substring(0,ReadSetArray_right[z].length()-((int)(0.2*ReadSetArray_right[z].length())))+"\n");
					   writer4.close();
				   }
				   if(ReadSetArray_left[z].charAt(0)!='#' && ReadSetArray_right[z].charAt(0)=='#'){
					   FileWriter writer1= new FileWriter(ParentPath+"/ReadFiles/LowDepthReads_left.fa",true);
					   writer1.write(">Node_"+(IndexLowDepthRead_left++)+"\n"+ReadSetArray_left[z]+"\n");
					   writer1.close();
					   FileWriter writer2= new FileWriter(ParentPath +"/ReadFiles/LowDepthReads_right.fa",true);
					   writer2.write(">Node_"+(IndexLowDepthRead_right++)+"\n"+ReadSetArray_right[z].substring(1,ReadSetArray_right[z].length())+"\n");
					   writer2.close();
					   FileWriter writer3= new FileWriter(ParentPath+"/ReadFiles/NormalDepthReads_left.fa",true);
					   writer3.write(">Node_"+(IndexNormalDepthRead_left++)+"\n"+ReadSetArray_left[z].substring(0,ReadSetArray_left[z].length()-((int)(0.2*ReadSetArray_left[z].length())))+"\n");
					   writer3.close();
					   FileWriter writer4= new FileWriter(ParentPath +"/ReadFiles/NormalDepthReads_right.fa",true);
					   writer4.write(">Node_"+(IndexNormalDepthRead_right++)+"\n"+ReplaceString.substring(0,ReplaceString.length()-((int)(0.2*ReplaceString.length())))+"\n");
					   writer4.close();
				   }
			   }
		  }
	      //Merge FASTA.      
		  Process p_MergeLow=null;
		  Runtime r_MergeLow=Runtime.getRuntime();
		  Process p_MergeNormal=null;
		  Runtime r_MergeNormal=Runtime.getRuntime();
	      try{
		  	    String[] cmd_MergeLow = { "sh", "-c", ParentPath+"/tools/shuffleSequences_fasta.pl "+ParentPath+"/ReadFiles/LowDepthReads_left.fa"+" "+ParentPath +"/ReadFiles/LowDepthReads_right.fa"+" "+ParentPath+"/ReadFiles/LowDepthReads.fa"};
		  	    p_MergeLow=r_MergeLow.exec(cmd_MergeLow);
				p_MergeLow.waitFor();	
		  	    String[] cmd_MergeNormal = { "sh", "-c", ParentPath+"/tools/shuffleSequences_fasta.pl "+ParentPath+"/ReadFiles/NormalDepthReads_left.fa"+" "+ParentPath +"/ReadFiles/NormalDepthReads_right.fa"+" "+ParentPath+"/ReadFiles/NormalDepthReads.fa"};
		  	    p_MergeNormal=r_MergeNormal.exec(cmd_MergeNormal);
				p_MergeNormal.waitFor();
	      }
	      catch(Exception e){
	    	    System.out.println("Step08 Error:"+e.getMessage());
	    	    e.printStackTrace();
	      }
		  //Free Array.
		  LowFreKmerArray=null;
		  LowFreHashTable_A=null;
		  LowFreHashTable_T=null;
		  LowFreHashTable_G=null;
		  LowFreHashTable_C=null;
		  ReadSetArray_left=null;
		  ReadSetArray_right=null;
		  long orz_rc = Math.abs(startMem_rc - r_rc.freeMemory());
		  long endTime_rc = System.currentTimeMillis();
		  System.out.println(" Time consumption:"+(endTime_rc-startTime_rc)+"ms. Memory consumption:"+orz_rc/1000000000+"GB] ");   
	      //Step09-Generating multiple sets of paired-end reads of different insert size.
	      System.out.print("Step09: Generating multiple sets of Paired-end reads");
	      long startTime_gp = System.currentTimeMillis(); 
	      Runtime r_gp = Runtime.getRuntime();
	      long startMem_gp = r_gp.freeMemory();
          //Definition of different insert size.
	      int insertsize1=150;
	      int sd_insertsize1=15;
	      int insertsize2=300;
	      int sd_insertsize2=30;
	      int insertsize3=500;
	      int sd_insertsize3=50;
	      int insertsize4=700;
	      int sd_insertsize4=70;
	      int insertsize5=900;
	      int sd_insertsize5=90;
	      int insertsize6=1500;
	      int sd_insertsize6=150;
	      int insertsize7=3000;
	      int sd_insertsize7=300;
	      int insertsize8=5000;
	      int sd_insertsize8=500;
	      int insertsize9=7500;
	      int sd_insertsize9=750;
	      int insertsize10=9500;
	      int sd_insertsize10=950;
		  Process p_art_1=null;
		  Runtime r_art_1=Runtime.getRuntime();
		  Process p_art_2=null;
		  Runtime r_art_2=Runtime.getRuntime();
		  Process p_art_3=null;
		  Runtime r_art_3=Runtime.getRuntime();
		  Process p_art_4=null;
		  Runtime r_art_4=Runtime.getRuntime();
		  Process p_art_5=null;
		  Runtime r_art_5=Runtime.getRuntime();
		  Process p_art_6=null;
		  Runtime r_art_6=Runtime.getRuntime();
		  Process p_art_7=null;
		  Runtime r_art_7=Runtime.getRuntime();
		  Process p_art_8=null;
		  Runtime r_art_8=Runtime.getRuntime();
		  Process p_art_9=null;
		  Runtime r_art_9=Runtime.getRuntime();
		  Process p_art_10=null;
		  Runtime r_art_10=Runtime.getRuntime();
	      try{
	    	    String[] cmd_artRun_1 = { "sh", "-c", ParentPath+"/tools/art_bin_MountRainier/art_illumina -ss HS25 -sam -i "+ParentPath+"/Assembly/SPAdes/scaffolds.ChangeLines.fasta -p -ef  -l 49 -f 100 -m "+insertsize1+" -s "+sd_insertsize1+" -o "+ParentPath+"/Auxiliary_lib/InsertSize"+insertsize1+"_  > "+ParentPath+"/Log/art.log"};
	    	    String[] cmd_artRun_2 = { "sh", "-c", ParentPath+"/tools/art_bin_MountRainier/art_illumina -ss HS25 -sam -i "+ParentPath+"/Assembly/SPAdes/scaffolds.ChangeLines.fasta -p -ef  -l 49 -f 100 -m "+insertsize2+" -s "+sd_insertsize2+" -o "+ParentPath+"/Auxiliary_lib/InsertSize"+insertsize2+"_  > "+ParentPath+"/Log/art.log"};
	    	    String[] cmd_artRun_3 = { "sh", "-c", ParentPath+"/tools/art_bin_MountRainier/art_illumina -ss HS25 -sam -i "+ParentPath+"/Assembly/SPAdes/scaffolds.ChangeLines.fasta -p -ef  -l 49 -f 100 -m "+insertsize3+" -s "+sd_insertsize3+" -o "+ParentPath+"/Auxiliary_lib/InsertSize"+insertsize3+"_  > "+ParentPath+"/Log/art.log"};
	    	    String[] cmd_artRun_4 = { "sh", "-c", ParentPath+"/tools/art_bin_MountRainier/art_illumina -ss HS25 -sam -i "+ParentPath+"/Assembly/SPAdes/scaffolds.ChangeLines.fasta -p -ef  -l 49 -f 100 -m "+insertsize4+" -s "+sd_insertsize4+" -o "+ParentPath+"/Auxiliary_lib/InsertSize"+insertsize4+"_  > "+ParentPath+"/Log/art.log"};
	    	    String[] cmd_artRun_5 = { "sh", "-c", ParentPath+"/tools/art_bin_MountRainier/art_illumina -ss HS25 -sam -i "+ParentPath+"/Assembly/SPAdes/scaffolds.ChangeLines.fasta -p -ef  -l 49 -f 100 -m "+insertsize5+" -s "+sd_insertsize5+" -o "+ParentPath+"/Auxiliary_lib/InsertSize"+insertsize5+"_  > "+ParentPath+"/Log/art.log"};
	    	    String[] cmd_artRun_6 = { "sh", "-c", ParentPath+"/tools/art_bin_MountRainier/art_illumina -ss HS25 -sam -i "+ParentPath+"/Assembly/SPAdes/scaffolds.ChangeLines.fasta -p -ef  -l 49 -f 100 -m "+insertsize6+" -s "+sd_insertsize6+" -o "+ParentPath+"/Auxiliary_lib/InsertSize"+insertsize6+"_  > "+ParentPath+"/Log/art.log"};
	    	    String[] cmd_artRun_7 = { "sh", "-c", ParentPath+"/tools/art_bin_MountRainier/art_illumina -ss HS25 -sam -i "+ParentPath+"/Assembly/SPAdes/scaffolds.ChangeLines.fasta -p -ef  -l 49 -f 100 -m "+insertsize7+" -s "+sd_insertsize7+" -o "+ParentPath+"/Auxiliary_lib/InsertSize"+insertsize7+"_  > "+ParentPath+"/Log/art.log"};
	    	    String[] cmd_artRun_8 = { "sh", "-c", ParentPath+"/tools/art_bin_MountRainier/art_illumina -ss HS25 -sam -i "+ParentPath+"/Assembly/SPAdes/scaffolds.ChangeLines.fasta -p -ef  -l 49 -f 100 -m "+insertsize8+" -s "+sd_insertsize8+" -o "+ParentPath+"/Auxiliary_lib/InsertSize"+insertsize8+"_  > "+ParentPath+"/Log/art.log"};
	    	    String[] cmd_artRun_9 = { "sh", "-c", ParentPath+"/tools/art_bin_MountRainier/art_illumina -ss HS25 -sam -i "+ParentPath+"/Assembly/SPAdes/scaffolds.ChangeLines.fasta -p -ef  -l 49 -f 100 -m "+insertsize9+" -s "+sd_insertsize9+" -o "+ParentPath+"/Auxiliary_lib/InsertSize"+insertsize9+"_  > "+ParentPath+"/Log/art.log"};
	    	    String[] cmd_artRun_10 = { "sh", "-c", ParentPath+"/tools/art_bin_MountRainier/art_illumina -ss HS25 -sam -i "+ParentPath+"/Assembly/SPAdes/scaffolds.ChangeLines.fasta -p -ef -l 49 -f 100 -m "+insertsize10+" -s "+sd_insertsize10+" -o "+ParentPath+"/Auxiliary_lib/InsertSize"+insertsize10+"_  > "+ParentPath+"/Log/art.log"};
	    	    p_art_1=r_art_1.exec(cmd_artRun_1);
	    	    p_art_2=r_art_2.exec(cmd_artRun_2);	
	    	    p_art_3=r_art_3.exec(cmd_artRun_3);
	    	    p_art_4=r_art_4.exec(cmd_artRun_4);
	    	    p_art_5=r_art_5.exec(cmd_artRun_5);
	    	    p_art_6=r_art_6.exec(cmd_artRun_6);
	    	    p_art_7=r_art_7.exec(cmd_artRun_7);
	    	    p_art_8=r_art_8.exec(cmd_artRun_8);
	    	    p_art_9=r_art_9.exec(cmd_artRun_9);
	    	    p_art_10=r_art_10.exec(cmd_artRun_10);
				p_art_1.waitFor();
				p_art_2.waitFor();
				p_art_3.waitFor();
				p_art_4.waitFor();
				p_art_5.waitFor();
				p_art_6.waitFor();
				p_art_7.waitFor();
				p_art_8.waitFor();
				p_art_9.waitFor();
				p_art_10.waitFor();
	      }
	      catch(Exception e){
	    	    System.out.println("Step09 Error:"+e.getMessage());
	    	    e.printStackTrace();
	      }
	      CommonClass.GenerateFastaFromFastqFiles(ParentPath,ParentPath+"/Auxiliary_lib/InsertSize"+insertsize1+"_1.fq",ParentPath+"/Auxiliary_lib/InsertSize"+insertsize1+"_2.fq","InsertSize"+insertsize1);
	      CommonClass.GenerateFastaFromFastqFiles(ParentPath,ParentPath+"/Auxiliary_lib/InsertSize"+insertsize2+"_1.fq",ParentPath+"/Auxiliary_lib/InsertSize"+insertsize2+"_2.fq","InsertSize"+insertsize2);
	      CommonClass.GenerateFastaFromFastqFiles(ParentPath,ParentPath+"/Auxiliary_lib/InsertSize"+insertsize3+"_1.fq",ParentPath+"/Auxiliary_lib/InsertSize"+insertsize3+"_2.fq","InsertSize"+insertsize3);
	      CommonClass.GenerateFastaFromFastqFiles(ParentPath,ParentPath+"/Auxiliary_lib/InsertSize"+insertsize4+"_1.fq",ParentPath+"/Auxiliary_lib/InsertSize"+insertsize4+"_2.fq","InsertSize"+insertsize4);
	      CommonClass.GenerateFastaFromFastqFiles(ParentPath,ParentPath+"/Auxiliary_lib/InsertSize"+insertsize5+"_1.fq",ParentPath+"/Auxiliary_lib/InsertSize"+insertsize5+"_2.fq","InsertSize"+insertsize5);
	      CommonClass.GenerateFastaFromFastqFiles(ParentPath,ParentPath+"/Auxiliary_lib/InsertSize"+insertsize6+"_1.fq",ParentPath+"/Auxiliary_lib/InsertSize"+insertsize6+"_2.fq","InsertSize"+insertsize6);
	      CommonClass.GenerateFastaFromFastqFiles(ParentPath,ParentPath+"/Auxiliary_lib/InsertSize"+insertsize7+"_1.fq",ParentPath+"/Auxiliary_lib/InsertSize"+insertsize7+"_2.fq","InsertSize"+insertsize7);
	      CommonClass.GenerateFastaFromFastqFiles(ParentPath,ParentPath+"/Auxiliary_lib/InsertSize"+insertsize8+"_1.fq",ParentPath+"/Auxiliary_lib/InsertSize"+insertsize8+"_2.fq","InsertSize"+insertsize8);
	      CommonClass.GenerateFastaFromFastqFiles(ParentPath,ParentPath+"/Auxiliary_lib/InsertSize"+insertsize9+"_1.fq",ParentPath+"/Auxiliary_lib/InsertSize"+insertsize9+"_2.fq","InsertSize"+insertsize9);
	      CommonClass.GenerateFastaFromFastqFiles(ParentPath,ParentPath+"/Auxiliary_lib/InsertSize"+insertsize10+"_1.fq",ParentPath+"/Auxiliary_lib/InsertSize"+insertsize10+"_2.fq","InsertSize"+insertsize10);
	      long orz_gp = Math.abs(startMem_gp - r_gp.freeMemory());
	      long endTime_gp = System.currentTimeMillis();
	      System.out.println(" [Time consumption:"+(endTime_gp-startTime_gp)+"ms. Memory consumption:"+orz_gp/1000000000+"GB] ");
	      //Step10-Assembling low depth Paired-end reads.
	      System.out.print("Step10: Assembling low depth Paired-end reads");
	      long startTime_AssemblyLow = System.currentTimeMillis(); 
	      Runtime r_AssemblyLow = Runtime.getRuntime();
	      long startMem_AssemblyLow = r_AssemblyLow.freeMemory();
	      //Assembly start.
		  Process p_assemblylow=null;
		  Runtime r_assemblylow=Runtime.getRuntime();
	      try{
	    	    String[] cmd_assemblylow = { "sh", "-c", ParentPath+"/tools/epga "+ParentPath+"/ReadFiles/LowDepthReads.fa "+i+" "+s+" "+k+" 16"};
	    	    p_assemblylow=r_assemblylow.exec(cmd_assemblylow);
				p_assemblylow.waitFor();
	      }
	      catch(Exception e){
	    	    System.out.println("Step10 Error:"+e.getMessage());
	    	    e.printStackTrace();
	      }
		  File LowDepthContigFiles=new File(ParentPath+"/contigSetLong.fa");
		  File LowDepthScaffoldFiles=new File(ParentPath+"/scaffoldLong.fa");
	      if(LowDepthContigFiles.exists())
	      {
	    	  CommonClass.copyFile(ParentPath+"/contigSetLong.fa", ParentPath+"/Assembly/EPGA-SC/LowDepthReads/contigSetLong.fa");
	      }
	      if(LowDepthScaffoldFiles.exists())
	      {
	    	  CommonClass.copyFile(ParentPath+"/scaffoldLong.fa", ParentPath+"/Assembly/EPGA-SC/LowDepthReads/scaffoldLong.fa");
	      }
	      CommonClass.delSpecialFile(directory,"DATA",".fa");
		  //Assembly end.
	      long orz_AssemblyLow = Math.abs(startMem_AssemblyLow - r_AssemblyLow.freeMemory());
	      long endTime_AssemblyLow = System.currentTimeMillis();
	      System.out.println(" [Time consumption:"+(endTime_AssemblyLow-startTime_AssemblyLow)+"ms. Memory consumption:"+orz_AssemblyLow/1000000000+"GB] ");
	      //Step11-Assembling normal depth Paired-end reads .
	      System.out.print("Step11: Assembling normal depth Paired-end reads");
		  long startTime_AssemblyNormal = System.currentTimeMillis(); 
	      Runtime r_AssemblyNormal = Runtime.getRuntime();
	      long startMem_AssemblyNormal = r_AssemblyNormal.freeMemory();
	      //Assembly start.
		  Process p_assemblynormal=null;
		  Runtime r_assemblynormal=Runtime.getRuntime();
	      try{
	    	    String[] cmd_assemblynormal = { "sh", "-c", ParentPath+"/tools/epga "+ParentPath+"/ReadFiles/NormalDepthReads.fa "+i+" "+s+" "+ParentPath+"/ReadFiles/InsertSize"+insertsize1+".fa "+insertsize1+" "+sd_insertsize1+" "+ParentPath+"/ReadFiles/InsertSize"+insertsize2+".fa "+insertsize2+" "+sd_insertsize2+" "+ParentPath+"/ReadFiles/InsertSize"+insertsize3+".fa "+insertsize3+" "+sd_insertsize3+" "+ParentPath+"/ReadFiles/InsertSize"+insertsize4+".fa "+insertsize4+" "+sd_insertsize4+" "+ParentPath+"/ReadFiles/InsertSize"+insertsize5+".fa "+insertsize5+" "+sd_insertsize5+" "+ParentPath+"/ReadFiles/InsertSize"+insertsize6+".fa "+insertsize6+" "+sd_insertsize6+" "+ParentPath+"/ReadFiles/InsertSize"+insertsize7+".fa "+insertsize7+" "+sd_insertsize7+" "+ParentPath+"/ReadFiles/InsertSize"+insertsize8+".fa "+insertsize8+" "+sd_insertsize8+" "+ParentPath+"/ReadFiles/InsertSize"+insertsize9+".fa "+insertsize9+" "+sd_insertsize9+" "+ParentPath+"/ReadFiles/InsertSize"+insertsize10+".fa "+insertsize10+" "+sd_insertsize10+" "+K+" 16"};
	    	    p_assemblynormal=r_assemblynormal.exec(cmd_assemblynormal);
				p_assemblynormal.waitFor();
	      }
	      catch(Exception e){
	    	    System.out.println("Step11 Error:"+e.getMessage());
	    	    e.printStackTrace();
	      }
		  File NormalDepthContigFiles=new File(ParentPath+"/contigSetLong.fa");
		  File NormalDepthScaffoldFiles=new File(ParentPath+"/scaffoldLong.fa");
	      if(NormalDepthContigFiles.exists())
	      {
	    	  CommonClass.copyFile(ParentPath+"/contigSetLong.fa", ParentPath+"/Assembly/EPGA-SC/NormalDepthReads/contigSetLong.fa");
	      }
	      if(NormalDepthScaffoldFiles.exists())
	      {
	    	  CommonClass.copyFile(ParentPath+"/scaffoldLong.fa", ParentPath+"/Assembly/EPGA-SC/NormalDepthReads/scaffoldLong.fa");
	      }
	      CommonClass.delSpecialFile(directory,"DATA",".fa");
		  //Assembly end.
		  long orz_AssemblyNormal = Math.abs(startMem_AssemblyNormal - r_AssemblyNormal.freeMemory());
	      long endTime_AssemblyNormal = System.currentTimeMillis();
	      System.out.println(" [Time consumption:"+(endTime_AssemblyNormal-startTime_AssemblyNormal)+"ms. Memory consumption:"+(double)orz_AssemblyNormal/1000000000+"GB] ");
	      System.out.print("Step12: Removing chimeric errors in assemblies");
	      long startTime_chimeric = System.currentTimeMillis(); 
	      Runtime r_chimeric = Runtime.getRuntime();
	      long startMem_chimeric = r_chimeric.freeMemory();
		  //Process.
	      Process pk_index=null;
	      Process pk_mapping=null;
	      Process pk_sam2bam=null;
		  //Runtime.
	      Runtime k_index = Runtime.getRuntime();
	      Runtime k_mapping = Runtime.getRuntime();
	      Runtime k_sam2bam  = Runtime.getRuntime();
	      //gun-ZIP.
		  Process p_gzip1=null;
		  Process p_gzip2=null;
		  Runtime r_gzip1=Runtime.getRuntime();
		  Runtime r_gzip2=Runtime.getRuntime();
	      try{
	    	    String[] cmd_gzip1 = { "sh", "-c", "gunzip "+ParentPath+"/Assembly/SPAdes/corrected/"+DataName1+".00.0_0.cor.fastq.gz"};
	    	    p_gzip1=r_gzip1.exec(cmd_gzip1);
				p_gzip1.waitFor();
				String[] cmd_gzip2 = { "sh", "-c", "gunzip "+ParentPath+"/Assembly/SPAdes/corrected/"+DataName2+".00.0_0.cor.fastq.gz"};
				p_gzip2=r_gzip2.exec(cmd_gzip2);
				p_gzip2.waitFor();
	      }
	      catch(Exception e){
	    	    System.out.println("Step14 Scaffolding pre-process:"+e.getMessage());
	    	    e.printStackTrace();
	      }
	      //Remove MUMmer files.
	      File SamFiles=new File(ParentPath+"/Alignments/");
	      File BamFiles=new File(ParentPath+"/Alignments/");
	      CommonClass.delSpecialFile(SamFiles,"Data_nucmer","DATA");
	      CommonClass.delSpecialFile(SamFiles,"Data_filter","DATA");
	      CommonClass.delSpecialFile(SamFiles,"DATA",".delta");
	      CommonClass.delSpecialFile(BamFiles,"DATA",".txt");
	      CommonClass.delSpecialFile(BamFiles,"MUMmerGroupRecords",".fa");
	      try{
	    	    //Index.
	    	    String[] cmd_nucmer = { "sh", "-c", ParentPath+"/tools/MUMmer3.23/nucmer -c 50 -p "+ParentPath+"/Alignments/Data_nucmer "+ParentPath+"/Assembly/EPGA-SC/NormalDepthReads/scaffoldLong.fa "+ParentPath+"/Assembly/SPAdes/scaffolds.ChangeLines.fasta"};
	    	    pk_index=k_index.exec(cmd_nucmer);
	    	    pk_index.waitFor();
				//Mapping.
	    	    String[] cmd_mapping = { "sh", "-c", ParentPath+"/tools/MUMmer3.23/delta-filter -i 97 -q -r "+ParentPath+"/Alignments/Data_nucmer.delta > "+ParentPath+"/Alignments/Data_filter"};
	    	    pk_mapping=k_mapping.exec(cmd_mapping);
				pk_mapping.waitFor();
				//SAM TO BAM.
	    	    String[] cmd_sam2bam = { "sh", "-c", ParentPath+"/tools/MUMmer3.23/show-coords -dTlro "+ParentPath+"/Alignments/Data_filter > "+ParentPath+"/Alignments/Data_coords.txt"};
	    	    pk_sam2bam=k_sam2bam.exec(cmd_sam2bam);
	    	    pk_sam2bam.waitFor();
	      }
	      catch(Exception e){
	    	    System.out.println("Step12 Align Error:"+e.getMessage());
	    	    e.printStackTrace();
	      }
	      //Break the miss-assembly.
  	      File Corr_contigs=new File(ParentPath+"/Assembly/EPGA-SC/FinalAssembly/");
  	      CommonClass.delSpecialFile(Corr_contigs,"contigs_corr",".fa");
	      CommonClass.breakErrorPoints(R,ParentPath,ParentPath+"/Alignments/Data_coords.txt",ParentPath+"/Assembly/EPGA-SC/NormalDepthReads/scaffoldLong.fa",ParentPath+"/Assembly/SPAdes/scaffolds.ChangeLines.fasta",ParentPath+"/Assembly/EPGA-SC/FinalAssembly","200");
	      long orz_chimeric = Math.abs(startMem_chimeric - r_chimeric.freeMemory());
	      long endTime_chimeric = System.currentTimeMillis();
	      System.out.println(" Time consumption:"+(endTime_chimeric-startTime_chimeric)+"ms. Memory consumption:"+(double)orz_chimeric/1000000000+"GB] ");
	      System.out.print("Step12: Merging and exstension of Assemblies");
	      long startTime_merging = System.currentTimeMillis(); 
	      Runtime r_merging = Runtime.getRuntime();
	      long startMem_merging = r_merging.freeMemory();
	      //Loading low depth assemblies.
	      int NumLow=0;
	      File LowDepthAssemblies=new File(ParentPath+"/Assembly/EPGA-SC/LowDepthReads/scaffoldLong.fa");
	      File NormalDepthAssemblies=new File(ParentPath+"/Assembly/EPGA-SC/FinalAssembly/contigs_corr.fa");
	      if(LowDepthAssemblies.exists()&&NormalDepthAssemblies.exists()){
		      int ArraySize_lowdepth=CommonClass.getFileLines(ParentPath+"/Assembly/EPGA-SC/LowDepthReads/scaffoldLong.fa")/2;
			  String[] ReadSetArray_lowdepth=new String[ArraySize_lowdepth];
			  int scount_lowdepth=CommonClass.FileToArray(ParentPath+"/Assembly/EPGA-SC/LowDepthReads/scaffoldLong.fa",ReadSetArray_lowdepth,">");			  
			  int ArraySize_normaldepth=CommonClass.getFileLines(ParentPath+"/Assembly/EPGA-SC/FinalAssembly/contigs_corr.fa")/2;
			  String[] ReadSetArray_normaldepth=new String[ArraySize_normaldepth];
			  int scount_normaldepth=CommonClass.FileToArray(ParentPath+"/Assembly/EPGA-SC/FinalAssembly/contigs_corr.fa",ReadSetArray_normaldepth,">");
			  if(scount_lowdepth>0&&scount_normaldepth>0){
				  for(int z=0;z<scount_lowdepth;z++){
					  int Flag=1;
			    	  for(int b=0;b<scount_normaldepth;b++){
			    		  if((ReadSetArray_normaldepth[b].contains(ReadSetArray_lowdepth[z]))||(ReadSetArray_normaldepth[b].contains(CommonClass.reverse(ReadSetArray_lowdepth[z])))){
			    			  Flag=0;
			    			  break;
			    		  }  
			    	  }
			    	  if(Flag==1){
		   	              FileWriter writer2= new FileWriter(ParentPath+"/Assembly/EPGA-SC/FinalAssembly/contigs_corr.fa",true);
		                  writer2.write(">L"+(NumLow++)+"\n"+ReadSetArray_lowdepth[z]+"\n");
		                  writer2.close();
			    	  }
			      }
			  }
		      ReadSetArray_lowdepth=null;
		      ReadSetArray_normaldepth=null;
	      }
	      long orz_merging = Math.abs(startMem_merging - r_merging.freeMemory());
	      long endTime_merging = System.currentTimeMillis();
	      System.out.println(" [Time consumption:"+(endTime_merging-startTime_merging)+"ms. Memory consumption:"+(double)orz_merging/1000000000+"GB] ");
	      System.out.print("Step14: Scaffolding");
	      long startTime_Scaffolding = System.currentTimeMillis(); 
	      Runtime r_Scaffolding = Runtime.getRuntime();
	      long startMem_Scaffolding = r_Scaffolding.freeMemory(); 
		  Process p_scaffolds=null;
		  Runtime r_scaffolds=Runtime.getRuntime();
		  //Write configure.
  	      File ScaffConfigFiles=new File(ParentPath+"/tools/SSPACE-STANDARD-3.0_linux-x86_64/");
  	      CommonClass.delSpecialFile(ScaffConfigFiles,"config",".txt");
		  String Command1="lib1\tbwa\t"+ParentPath+"/Assembly/SPAdes/corrected/"+DataName1+".00.0_0.cor.fastq "+ParentPath+"/Assembly/SPAdes/corrected/"+DataName2+".00.0_0.cor.fastq "+i+"\t"+0.1+"\t"+"FR"+"\n"+"lib2\tbwa\t"+ParentPath+"/ReadFiles/read_left.fa "+ParentPath+"/ReadFiles/read_right.fa "+i+"\t"+0.1+"\t"+"FR";
		  FileWriter writer= new FileWriter(ParentPath+"/tools/SSPACE-STANDARD-3.0_linux-x86_64/config.txt",true);
		  writer.write(Command1+"\n");
	      writer.close();
	      try{
	    	    String[] cmd_scaffolds={ "sh", "-c", "perl "+ParentPath+"/tools/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE_Standard_v3.0.pl -a 0.85 -l "+ParentPath+"/tools/SSPACE-STANDARD-3.0_linux-x86_64/config.txt"+" -s "+ParentPath+"/Assembly/EPGA-SC/FinalAssembly/contigs_corr.fa -b ./scaffolds_sspace -T "+t};
	    	    p_scaffolds=r_scaffolds.exec(cmd_scaffolds);
				p_scaffolds.waitFor();
	      }
	      catch(Exception e){
	    	    System.out.println("Step14 Scaffolding:"+e.getMessage());
	    	    e.printStackTrace();
	      }
	      //Output final assemblies.
		  File FinalContigFile=new File(ParentPath+"/Assembly/EPGA-SC/FinalAssembly/contigs_corr.fa");
		  File FinalScaffoldFile=new File(ParentPath+"/scaffolds_sspace/scaffolds_sspace.final.scaffolds.fasta");
	      if(FinalContigFile.exists())
	      {
	    	  CommonClass.copyFile(ParentPath+"/Assembly/EPGA-SC/FinalAssembly/contigs_corr.fa", o+"/Contigs.fa");
	      }
	      if(FinalScaffoldFile.exists())
	      {
	    	  CommonClass.copyFile(ParentPath+"/scaffolds_sspace/scaffolds_sspace.final.scaffolds.fasta", o+"/Scaffolds.fa");
	      }
	      long orz_Scaffolding = Math.abs(startMem_Scaffolding - r_Scaffolding.freeMemory());
	      long endTime_Scaffolding = System.currentTimeMillis();
	      System.out.println(" [Time consumption:"+(endTime_Scaffolding-startTime_Scaffolding)+"ms. Memory consumption:"+(double)orz_Scaffolding/1000000000+"GB] ");
	  }
	  else
	  {
		  System.out.println("\n -q1, -q2 and -o are three required parameters. \nPlease check the configuration of these three parameters.\n");
	  }
   }
}
