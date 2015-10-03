package com.antigenomics.vdjstruct;

import java.util.Vector;
import java.util.Collections;
import java.io.BufferedWriter;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.FileWriter;
import java.io.FileReader;
import com.milaboratory.core.sequence.AminoAcidSequence;

public class Table
{
	public int AlphaNum;
	public int BetaNum;
	
	public Vector<Chain> AlphaArray;
	public Vector<Chain> BetaArray;

	public Graph AlphaGraph;
	public Graph BetaGraph;
	
//===================================================================================
//================= Routine functions ===============================================
//===================================================================================

	public void ReadSampleTable() throws IOException
	{
		BufferedReader file = new BufferedReader(new FileReader("data/vdjdb.txt"));
		String line;
		file.readLine();
		int i = 0;
		int bflag = 0;
		int aflag = 0;
		BetaNum = 0;
		AlphaNum = 0;
		while((line = file.readLine()) != null)
		{
			String[] list = line.split("\t");
			if (list[3].equals("TRB"))
			{
				if (bflag == 0)
				{
					BetaArray = new Vector<>();
					bflag = 1;
				}
				if (!list[8].equals("."))
				{
					BetaArray.addElement(new Chain("TRB",
							".",
							i,
							".",
							new int[]{-2, -2},
							list[1],
							new int[]{-2, -2},
							list[2],
							new int[]{-2, -2},
							".",
							new int[]{-2, -2},
							".",
							new AminoAcidSequence(list[0].replaceAll(" ", "")),
							new AminoAcidSequence(list[8].replaceAll(" ", ""))));
					BetaNum++;
					i++;
				}
			}
			else
				if (list[3].equals("TRA"))
				{
					if (aflag == 0)
					{
						AlphaArray = new Vector<>();
						aflag = 1;
					}
					if (!list[8].equals("."))
					{
						AlphaArray.addElement(new Chain("TRA",
								".",
								i,
								".",
								new int[]{-2, -2},
								list[1],
								new int[]{-2, -2},
								list[2],
								new int[]{-2, -2},
								".",
								new int[]{-2, -2},
								".",
								new AminoAcidSequence(list[0].replaceAll(" ", "")),
								new AminoAcidSequence(list[8].replaceAll(" ", ""))));
						AlphaNum++;
						i++;
					}
				}
		}
		if (AlphaNum != 0)
		{
			Collections.sort(AlphaArray, Chain.sortByCDR());
			for (int j = 0; j < AlphaArray.size()-1; j++)
				if (Chain.sortByCDR().compare(AlphaArray.get(j), AlphaArray.get(j+1)) == 0) {
					AlphaArray.remove(j--);
					AlphaNum--;
				}
			Collections.sort(AlphaArray, Chain.sortById());
		}
		if (BetaNum != 0)
		{
			Collections.sort(BetaArray, Chain.sortByCDR());
			/*for (int k = 0; k < BetaArray.size(); k++)
				System.out.println(BetaArray.get(k).getCDR3().toString());*/
			for (int j = 0; j < BetaArray.size()-1; j++)
				if (Chain.sortByCDR().compare(BetaArray.get(j), BetaArray.get(j+1)) == 0) {
					BetaArray.remove(j--);
					BetaNum--;
				}
			Collections.sort(BetaArray, Chain.sortById());
		}
	}
	public void AlphaRead(String fileName) throws IOException
	{
		BufferedReader file = new BufferedReader(new FileReader(fileName));
		String line;
		AlphaArray = new Vector<>();
		int i = 0;
		int flag = 1;
		int[] array = new int[8];
		
		while((line = file.readLine()) != null)
		{
			String[] list = line.split("/");
			int buf = 0;
			
			if (list.length != 16)
				flag = 0;
			else
				for (int n = 0; n < 12; n++)
				{	
					if ((n+1)%3 != 0)
						array[n-buf] = Integer.parseInt(list[3+n]);
					else
						buf++;
						
					if (array[n-buf] == -1)
						if ((n-buf != 2) && (n-buf != 3) && (n-buf != 6) && (n-buf != 7))
						{
							flag = 0;
							break;
						}
				}
		
			if (flag != 0)
			{
				AlphaArray.addElement(new Chain("TRA",
					list[1],
						i,
					list[2], 
					new int[]{array[0], array[1]},
					list[5],
					new int[]{array[4], array[5]},
					list[11],
					new int[]{array[6], array[7]},
					list[14],
					new int[]{-2, -2},
					".",
						new AminoAcidSequence(list[15]),
						null));
				i++;
			}
			flag = 1;
		}
		AlphaNum = i;
		file.close();
	}
	
	public void BetaRead(String fileName) throws IOException
	{
		BufferedReader file = new BufferedReader(new FileReader(fileName));
		String line;
		BetaArray = new Vector<>();
		int i = 0;
		int flag = 1;
		int[] array = new int[8];
		while((line = file.readLine()) != null)
		{
			String[] list = line.split("/");
			int buf = 0;
			
			if (list.length != 16)
				flag = 0;
			else
				for (int n = 0; n < 12; n++)
				{	
					if ((n+1)%3 != 0)
						array[n-buf] = Integer.parseInt(list[3+n]);
					else
						buf++;
						
					if (array[n-buf] == -1)
						if ((n-buf != 2) && (n-buf != 3) && (n-buf != 6) && (n-buf != 7))
						{
							flag = 0;
							break;
						}
				}
		
			if (flag != 0)
			{
				BetaArray.addElement(new Chain("TRB",
					list[1],
						i,
					list[2], 
					new int[]{array[0], array[1]},
					list[5],
					new int[]{array[2], array[3]},
					list[8],
					new int[]{array[4], array[5]},
					list[11],
					new int[]{array[6], array[7]},
					list[14],
						new AminoAcidSequence(list[15]),
						null));
				i++;
			}
			flag = 1;
		}
		BetaNum = i;
		file.close();
	}
	
	public void AlphaWrite(String fileName) throws IOException
	{
		BufferedWriter file = new BufferedWriter(new FileWriter(fileName));
		int i;
		Chain buf;
		for (i = 0; i < AlphaNum; i++)
		{
			buf = AlphaArray.get(i);
			buf.Write(file);
			file.write("\n");
		}
		file.close();
	}
	
	public void BetaWrite(String fileName) throws IOException
	{
		BufferedWriter file = new BufferedWriter(new FileWriter(fileName));
		int i;
		Chain buf;
		for (i = 0; i < BetaNum; i++)
		{
			buf = BetaArray.get(i);
			buf.Write(file);
			file.write("\n");
		}
		file.close();
	}	
}
	
