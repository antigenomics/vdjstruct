package com.antigenomics.vdjstruct;

import com.milaboratory.core.alignment.AffineGapAlignmentScoring;
import com.milaboratory.core.alignment.Aligner;
import com.milaboratory.core.alignment.BLASTMatrix;
import com.milaboratory.core.alignment.LinearGapAlignmentScoring;
import com.milaboratory.core.sequence.AminoAcidSequence;
import org.junit.Test;

public class GlobalAlignmentTest {
    @Test
    @SuppressWarnings("unchecked")
    public void test() {
        AminoAcidSequence cdr3a = new AminoAcidSequence("CAGSLSWGGFYNEQFF"),
                cdr3b = new AminoAcidSequence("CARSWAGGPRDEQYF");
        
        /*
         * Linear gap penalty, penalty for total number of gaps
         */

        System.out.println("Linear gap m=1,mm=-3,e=-5");

        System.out.println(
                Aligner.alignGlobalLinear(
                        new LinearGapAlignmentScoring<>(AminoAcidSequence.ALPHABET, 1, -3, -5), cdr3a, cdr3b
                ).getAlignmentHelper()
        );


        /*
         * Specifying amino acid substitution matrix (PAM250) instead of match/mismatch scores
         */

        System.out.println("\nLinear gap PAM250,e=-5");

        System.out.println(
                Aligner.alignGlobalLinear(
                        new LinearGapAlignmentScoring<>(AminoAcidSequence.ALPHABET, BLASTMatrix.BLOSUM45.getMatrix(), -5), cdr3a, cdr3b
                ).getAlignmentHelper()
        );

        System.out.println("\nAffine gap m=1,mm=-3,o=-5,e=-1");

        System.out.println(
                Aligner.alignGlobalAffine(
                        new AffineGapAlignmentScoring(AminoAcidSequence.ALPHABET, 1, -4, -5, -1), cdr3a, cdr3b
                ).getAlignmentHelper()
        );

        System.out.println("\nAffine gap PAM50,o=-5,e=-1");

        System.out.println(
                Aligner.alignGlobalAffine(
                        new AffineGapAlignmentScoring(AminoAcidSequence.ALPHABET, BLASTMatrix.BLOSUM45.getMatrix(), -5, -1), cdr3a, cdr3b
                ).getAlignmentHelper()
        );
    }
}
