package com.antigenomics.vdjstruct;

import com.milaboratory.core.sequence.AminoAcidSequence;

import java.util.Vector;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.FileWriter;

public class Graph
{
	int[][] Data;
	int Size;
	
	public void GraphFill(Vector<Chain> Array, int flag)
	{
		if (Array != null)
			if (!Array.isEmpty())
			{
				Size = Array.size();
				Data = new int[Size][Size];
				switch (flag)
				{
					case 0:
						this.GraphFill0(Array);
						break;
					case 1:
						this.GraphFill1(Array);
						break;
					case 2: 
						this.HeatMap(Array);
						break;
					case 3:
						this.UpgradedHeatMap(Array);
						break;
				}
			}
			else
			{
				System.out.println("Table is empty!");
				System.exit(0);
			}
		else
		{
				System.out.println("Null pointer received by Graph!");
				System.exit(0);
		}
	}
	
	private void GraphFill0(Vector<Chain> Array)
	{
		for (int i = 0; i < Size; i++)
		{
			for (int j = i; j < Size; j++)
			{
				if (Array.get(i).getAntigen().equals(Array.get(j).getAntigen()))
				{
					Data[i][j] = 0;
					Data[j][i] = 0;
				}
				else
				{
					Data[i][j] = 1;
					Data[j][i] = 1;
				}
			}
		}
	}
	
	private void GraphFill1(Vector<Chain> Array)
	{
		for (int i = 0; i < Size; i++)
		{
			for (int j = i; j < Size; j++)
			{
				String CDR3_i = Array.get(i).getCDR3().toString();
				String CDR3_j = Array.get(j).getCDR3().toString();
				double dist = (double)Utilities.EditorialDistance(CDR3_i, CDR3_j);
				double thres = Variable.EditDistThesh;
				if ((dist/(double)CDR3_i.length() < thres) || (dist/(double)CDR3_j.length() < thres))
				{
					Data[i][j] = 1;
					Data[j][i] = 1;
				}
				else
				{
					Data[i][j] = 0;
					Data[j][i] = 0;
				}
			}
		}
	}
	
	private void HeatMap(Vector<Chain> Array)
	{
		for (int i = 0; i < Size; i++)
		{
			for (int j = i; j < Size; j++)
			{
				String CDR3_i = Array.get(i).getCDR3().toString();
				String CDR3_j = Array.get(j).getCDR3().toString();
				int dist = Utilities.EditorialDistance(CDR3_i, CDR3_j);
				int maxlen = Math.max(CDR3_i.length(), CDR3_j.length());
				int minlen = Math.min(CDR3_i.length(), CDR3_j.length());
				Data[i][j] = (int)(1000.0*((double)dist/maxlen+(double)dist/minlen));
				Data[j][i] = Data[i][j];
			}
		}
	}

	private void UpgradedHeatMap(Vector<Chain> Array)
	{
		for (int i = 0; i < Size; i++)
		{
			for (int j = i; j < Size; j++)
			{
				int len_i = Array.get(i).getCDR3().size();
				int len_j = Array.get(i).getCDR3().size();
				AminoAcidSequence CDR3_i = Array.get(i).getCDR3();
				AminoAcidSequence CDR3_j = Array.get(j).getCDR3();
				//System.out.print(Utilities.DetailedAlignment(CDR3_i, CDR3_j)+" ");
				Data[i][j] = (int)((10.0 - Utilities.DetailedAlignment(CDR3_i, CDR3_j)/(double)(len_i+len_j))*1000.0);
				Data[j][i] = Data[i][j];
			}
		}
	}
	
	public void PrintPseudoGraph(String fileName) throws IOException
	{
		BufferedWriter file = new BufferedWriter(new FileWriter(fileName));
		System.out.println("Graph size = "+Size);
		for (int i = 0; i < Size; i++)
		{
			for (int j = 0; j < Size; j++)
				if (i > j)
					file.write(" ");
				else
					if (Data[i][j] == 1)
						file.write(9608);
					else
						if (Data[i][j] == 0)
							file.write(" ");
						else
						{
							System.out.println("Incorrect symbol in Data!");
							System.exit(0);
						}
			file.write("\n");
		}
		file.close();
	}
	
	public void PrintHeatMap(String fileName) throws IOException
	{
		BufferedWriter file = new BufferedWriter(new FileWriter(fileName));
		System.out.println("Graph size = "+Size);
		for (int i = 0; i < Size; i++)
		{
			for (int j = 0; j < Size; j++)
			{
				file.write(Integer.toString(Data[i][j]));
				if (j != Size-1)
					file.write(" ");
			}
			file.write("\n");
		}
		file.close();
	}
}
