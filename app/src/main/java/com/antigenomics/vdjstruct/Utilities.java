package com.antigenomics.vdjstruct;

import java.util.Vector;
import java.util.Arrays;
import java.util.Collections;
import java.io.BufferedWriter;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.FileWriter;
import java.io.FileReader;

public class Utilities
{	
	public static int EditorialDistance(String a, String b)
	{
		int len_a = a.length();
		int len_b = b.length();
		boolean flag = false;
		int [][] cols = new int [2][len_a+1];
		int index1 = 0;
		int index2 = 1;
		for (int i = 0; i <= len_a; i++)
			cols[0][i] = i;
		for (int i = 1; i <= len_b; i++)
		{
			if (flag)
			{
				index2 = 0;
				index1 = 1;
			}
			else 
			{
				index2 = 1;
				index1 = 0;
			}
				
			for (int j = 0; j <= len_a; j++)
			{
				if (j == 0)
					cols[index2][0] = i;
				else
					if (b.charAt(i-1) == a.charAt(j-1))
						cols[index2][j] = cols[index1][j-1];
					else
						cols[index2][j] = Min3(cols[index1][j-1], cols[index1][j], cols[index2][j-1])+1;
						
			}
			flag = !flag;
		}
		return cols[index2][len_a];
	}
	
	private static int[][] strAlignmentDetailed(String a, String b)
	{
		int len_a = a.length();
		int len_b = b.length();
		int [][] cols = new int [len_a+1][len_b+1];
		for (int i = 0; i <= len_a; i++)
			cols[i][0] = i;
		for (int j = 1; j <= len_b; j++)	
			for (int i = 0; i <= len_a; i++)
			{
				if (i == 0)
					cols[0][j] = j;
				else
					if (b.charAt(j-1) == a.charAt(i-1))
						cols[i][j] = cols[i-1][j-1];
					else
						cols[i][j] = Min3(cols[i-1][j-1], cols[i-1][j], cols[i][j-1])+1;
						
			}
			
		int i = len_a;
		int j = len_b;
		int[][] result = new int[2][len_a];
		Arrays.fill(result[1], 0);
		
		while(i > 0 || j > 0)
		{
			if (i > 0 && j > 0)
				if (cols[i-1][j] < cols[i][j-1])
				{
					result[0][--i] = -1;
					if (cols[i][j-1] < cols[i][j]+1)
						if (a.charAt(i) == b.charAt(j-1))
						{
							result[0][i] = 0;
							j--;
						}
						else
							if (cols[i][j-1] < cols[i][j])
							{
								result[0][i] = 2;
								j--;
							}
							/*else
							{
								
							}*/
				}
				else
				{
					//result[j-1] = 1;
					j--;
					if (cols[i-1][j] < cols[i][j]+1)
						if (a.charAt(i-1) == b.charAt(j))
							result[0][--i] = 0;
						else
							if (cols[i-1][j] < cols[i][j])
								result[0][--i] = 2;
							else 
								result[1][i-1]++;
				}
			else
				if (i == 0)
				{
					result[1][i-1]++;
				}
				else
					result[0][--i] = -1;
				/*	
				for (int k = 0; k < result[0].length; k++)	
					System.out.print(result[0][k]+" ");
				System.out.println();
				for (int k = 0; k < result[1].length; k++)	
					System.out.print(result[1][k]+" ");
				System.out.println();	
				System.out.println();*/
		}
		return result; // -1: deletion; 0: nothing; 1: insertion; 2: substitution
	}
	
	public static double strDistanceUpgraded(String a, String b, double ins, double del, double sub)
	{
		double dist = 0.0;
		int[][] ar = strAlignmentDetailed(a, b);
		for (int i = 0; i < ar[0].length; i++)
		{
			switch (ar[0][i])
			{
				case -1: dist += del; break;
				case  2: dist += sub; break;
			}
			if (ar[1][i] != 0)
				dist += (double)ar[1][i]*ins;
		}
		return (dist/(double)a.length()+dist/(double)b.length())*1000.0;
	}
	
	public static int Min3(int a, int b, int c)
	{
		return Math.min(c, Math.min(a, b));
	}
	
	public static void MistakeEmptyCheck(Vector array, String s)
	{
		if (array == null)
		{
			System.out.println(s+": Array is NULL!");
			System.exit(0);
		}
		if (array.isEmpty())
		{
			System.out.println(s+": Array is Empty!");
			System.exit(0);
		}
	}
	
	public static void main(String[] args)
	{
		String s1 = new String("qwertkjeekjfheeeelkjhlkjhlkjhpjjjjj");
		String s2 = new String("qwertkj;lkj;ljklkj;;lkj;ll");
		/*
		int len_a = s1.length();
		int len_b = s2.length();
		int [][] cols = new int [len_a+1][len_b+1];
		for (int i = 0; i <= len_a; i++)
			cols[i][0] = i;
		for (int j = 1; j <= len_b; j++)	
			for (int i = 0; i <= len_a; i++)
			{
				if (i == 0)
					cols[0][j] = j;
				else
					if (s2.charAt(j-1) == s1.charAt(i-1))
						cols[i][j] = cols[i-1][j-1];
					else
						cols[i][j] = Min3(cols[i-1][j-1], cols[i-1][j], cols[i][j-1])+1;
						
			}*/
			
		//System.out.println("last: "+cols[len_a][len_b]);	
		
		
		
		int [][] a = strAlignmentDetailed(s1, s2);
		for (int i = 0; i < a[0].length; i++)
		{
			System.out.print(a[0][i]+"\t");
			System.out.println();
		}
		
		/*for (int i = 0; i < s1.length(); i++)
		{
			for (int j = 0; j < s2.length(); j++)
				System.out.print(cols[i][j]+" ");
			System.out.println();
		}*/
		System.out.println("UpDist: "+strDistanceUpgraded(s1, s2, 1.0, 1.0, 1.0));
		System.out.println("OldDist: "+EditorialDistance(s1, s2));
	}
}
