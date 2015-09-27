package com.antigenomics.vdjstruct;

import java.util.Vector;

public class Claster
{
	Vector<Chain> Group;
	
	Claster()
	{
		Group = new Vector<Chain>();
	}
	
	public void AddChain(Chain c)
	{
		Chain d = new Chain();
		c.Copy(d);
		Group.addElement(d);
	}
	
	static Vector CompareClusters(Claster a, Claster b)
	{
		Vector<Chain> result = new Vector<Chain>();
		String Index_a;
		String Index_b;
		for (int i = 0; i < a.Group.size(); i++)
		{
			Index_a = a.Group.get(i).Id;
			for (int j = 0; j < b.Group.size(); j++)
			{
				Index_b = b.Group.get(i).Id;
				if (Index_a.equals(Index_b))
				{
					result.addElement(a.Group.get(i));
					a.Group.remove(i);
					b.Group.remove(j);
					break;
				}
			}
		}
		return result;
	}
}
