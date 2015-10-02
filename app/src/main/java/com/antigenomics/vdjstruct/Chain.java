package com.antigenomics.vdjstruct;

import java.io.IOException;
import java.io.BufferedWriter;
import java.util.Comparator;
import com.milaboratory.core.sequence.AminoAcidSequence;

public class Chain
{
	private String Id;
	private int Number;
	
	private String ChainFlag;// TRB/TRA

	private String Seq;
	private int[] VReg;
	private String VRegName;
	private int[] JReg;
	private String JRegName;
	private int[] CReg;	// Optional
	private String CRegName;// Optional
	private int[] DReg;	// Optional
	private String DRegName;// Optional

	private AminoAcidSequence CDR3;
	private AminoAcidSequence Antigen;	// Optional
	
	public Chain(String CFlag,
				 String i,
				 int Number,
				 String S,
				 int[] V,
				 String Vname,
				 int[] D,
				 String Dname,
				 int[] J,
				 String Jname,
				 int[] C,
				 String Cname,
				 AminoAcidSequence CDR3,
				 AminoAcidSequence Antigen) {
		this.Id = i;
		this.Number = Number;
		this.ChainFlag = CFlag;
		this.Seq = S;
		this.VReg = V;
		this.JReg = J;
		this.CReg = C;
		this.DReg = D;
		this.CDR3 = CDR3;
		this.VRegName = Vname;
		this.JRegName = Jname;
		this.CRegName = Cname;
		this.DRegName = Dname;
		this.Antigen = Antigen;
	}

	public Chain(Chain that)
	{
		if (that == null)
			throw new NullPointerException();

		this.ChainFlag 	= new String(that.ChainFlag);
		this.Id 		= new String(that.Id);
		this.Number 	= that.Number;
		this.Seq 		= new String(that.Seq);
		this.VReg 		= that.VReg.clone();
		this.VRegName 	= new String(that.VRegName);
		this.JReg 		= that.JReg.clone();
		this.JRegName 	= new String(that.JRegName);
		this.CReg 		= that.CReg.clone();
		this.CRegName 	= new String(that.CRegName);
		this.DReg 		= that.DReg.clone();
		this.DRegName 	= new String(that.DRegName);
		this.CDR3 		= new AminoAcidSequence(this.CDR3.asArray());
		this.Antigen 	= new AminoAcidSequence(this.Antigen.asArray());
	}

	public Chain Copy() {
		return new Chain(this);
	}
	
	public void Write(BufferedWriter file) throws IOException {
		file.write(this.toString());
	}
	
	public String toString() {
		return "Chain = " + ChainFlag +
				"\nId = " + Id +
				"\nSeq = " + Seq +
				"\nVReg = " + VReg[0] + "|" + VReg[1] +
				"\nVRegName = " + VRegName +
				"\nDReg = " + DReg[0] + "|" + DReg[1] +
				"\nDRegName = " + DRegName +
				"\nJReg = " + JReg[0] + "|" + JReg[1] +
				"\nJRegName = " + JRegName +
				"\nCReg = " + CReg[0] + "|" + CReg[1] +
				"\nCRegName = " + CRegName +
				"\nCDR3 = " + CDR3 +
				"\nAntigen = " + Antigen;
	}
	
	public static Comparator<Chain> sortByCDR()
	{
		return new ByCDR();
	}
	
	private static class ByCDR implements Comparator<Chain>
	{
		public int compare(Chain a, Chain b)
		{
			return a.CDR3.toString().compareTo(b.CDR3.toString());
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
			return a.Antigen.toString().compareTo(b.Antigen.toString());
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
	
	public String getD() {
		return Seq.substring(DReg[0], DReg[1]);
	}
	
	public String getSeq()
	{
		return Seq;
	}
	
	public String getId()
	{
		return Id;
	}

	public AminoAcidSequence getCDR3() {
		return CDR3;
	}

	public AminoAcidSequence getAntigen() {
		return Antigen;
	}
}
