package com.antigenomics.vdjstruct;

import java.io.IOException;
import java.io.BufferedWriter;
import java.util.Comparator;

public class Chain
{
	String Id;
	int Number;
	
	String ChainFlag;// TRB/TRA
	
	String Seq;
	int[] VReg;
	String VRegName;
	int[] JReg;
	String JRegName;
	int[] CReg;	// Optional
	String CRegName;// Optional
	int[] DReg;	// Optional
	String DRegName;// Optional
	
	String CDR3;
	String Antigen;	// Optional
	
	Chain()
	{
		ChainFlag 	= new String();
		Id 		= new String();
		Seq 		= new String();
		VReg 		= new int[2];
		VRegName 	= new String();
		JReg 		= new int[2];
		JRegName 	= new String();
		CReg 		= new int[2];
		CRegName 	= new String();
		DReg 		= new int[2];
		DRegName 	= new String();
		CDR3 		= new String();
		Antigen 	= new String();
	}
	
	Chain(String CFlag, String i, String S, int[] V, String Vname, int[] D,  String Dname, int[] J,  String Jname, int[] C,  String Cname, String CDR)
	{
		Id = i;
		ChainFlag = CFlag;
		Seq = S;
		VReg = V;
		JReg = J;
		CReg = C;
		DReg = D;
		CDR3 = CDR;
		VRegName = Vname;
		JRegName = Jname;
		CRegName = Cname;
		DRegName = Dname;
		Antigen = "unknown";
	}
	
	public void Write(BufferedWriter file) throws IOException
	{
		file.write("Chain = "+ChainFlag+
			"\nId = "+Id+
			"\nSeq = "+Seq+
			"\nVReg = "+VReg[0]+"|"+VReg[1]+
			"\nVRegName = "+VRegName+
			"\nDReg = "+DReg[0]+"|"+DReg[1]+
			"\nDRegName = "+DRegName+
			"\nJReg = "+JReg[0]+"|"+JReg[1]+
			"\nJRegName = "+JRegName+
			"\nCReg = "+CReg[0]+"|"+CReg[1]+
			"\nCRegName = "+CRegName+
			"\nCDR3 = "+CDR3+
			"\nAntigen = "+Antigen);
	}
	
	public void Show()
	{
		System.out.println("Chain = "+ChainFlag+
			"\nId = "+Id+
			"\nSeq = "+Seq+
			"\nVReg = "+VReg[0]+"|"+VReg[1]+
			"\nVRegName = "+VRegName+
			"\nDReg = "+DReg[0]+"|"+DReg[1]+
			"\nDRegName = "+DRegName+
			"\nJReg = "+JReg[0]+"|"+JReg[1]+
			"\nJRegName = "+JRegName+
			"\nCReg = "+CReg[0]+"|"+CReg[1]+
			"\nCRegName = "+CRegName+
			"\nCDR3 = "+CDR3+
			"\nAntigen = "+Antigen);
	}

	public void Copy(Chain target)
	{
		if (target != null)
		{
			target = new Chain();
		}
		target.ChainFlag 	= new String(this.ChainFlag);
		target.Id 		= new String(this.Id);
		target.Number 		= this.Number;
		target.Seq 		= new String(this.Seq);
		target.VReg 		= (int[])this.VReg.clone();
		target.VRegName 	= new String(this.VRegName);
		target.JReg 		= (int[])this.JReg.clone();
		target.JRegName 	= new String(this.JRegName);
		target.CReg 		= (int[])this.CReg.clone();
		target.CRegName 	= new String(this.CRegName);
		target.DReg 		= (int[])this.DReg.clone();
		target.DRegName 	= new String(this.DRegName);
		target.CDR3 		= new String(this.CDR3);
		target.Antigen 		= new String(this.Antigen);
	}
	
	public static Comparator<Chain> sortByCDR()
	{
		return new ByCDR();
	}
	
	private static class ByCDR implements Comparator<Chain>
	{
		public int compare(Chain a, Chain b)
		{
			return a.CDR3.compareTo(b.CDR3);
		}
	}
	
	public static Comparator<Chain> sortById()
	{
		return new ById();
	}
	
	private static class ById implements Comparator<Chain>
	{
		public int compare(Chain a, Chain b)
		{
			return a.Number - b.Number;
		}
	}
	
	public static Comparator<Chain> sortByAntigen()
	{
		return new ByAntigen();
	}
	
	private static class ByAntigen implements Comparator<Chain>
	{
		public int compare(Chain a, Chain b)
		{
			return a.Antigen.compareTo(b.Antigen);
		}
	}
	
	public String getV()
	{
		return Seq.substring(VReg[0], VReg[1]);
	}
	
	public String getJ()
	{
		return Seq.substring(JReg[0], JReg[1]);
	}
	
	public String getC()
	{
		return Seq.substring(CReg[0], CReg[1]);
	}
	
	public String getD()
	{
		return Seq.substring(DReg[0], DReg[1]);
	}
	
	public String getCDR3()
	{
		return CDR3;
	}
	
	public String getSeq()
	{
		return Seq;
	}
	
	public String getId()
	{
		return Id;
	}
}
