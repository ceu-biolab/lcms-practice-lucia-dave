package adduct;

import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/* Adduct: product of a direct addition of two or more distinct molecules,
 * resulting in a single reaction product containing all atoms of all components.
 */
public class Adduct {

    /**
     * Calculate the mass to search depending on the adduct hypothesis.
     * Supports single and multiple charges, as well as multimers (e.g., [2M+Na]+).
     *
     * @param mz mz
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return the monoisotopic mass of the experimental mass mz with the adduct @param adduct
     */
    public static Double getMonoisotopicMassFromMZ(Double mz, Map.Entry<String, Double> adduct) {

        Double monoisotopicMass = null;
        Double adductMass = adduct.getValue();
        int multimer = extractMultimer(adduct.getKey());
        int charge = extractCharge(adduct.getKey());

        // CASE 1: Single charge (no effect from charge or multimer)
        // If Adduct is single charge the formula is M = m/z +- adductMass. Charge is 1 so it does not affect
        if (charge == 1 && multimer == 1) {
            monoisotopicMass = mz + adductMass;
            //System.out.println(" -M = " + mz + " + " + adductMass + " = " + monoisotopicMass);
        }
        // CASE 2: Multiple charges, e.g., [M+2H]2+
        // If Adduct is double or triple charged the formula is M = ( mz +- adductMass ) * charge
        else if (charge > 1 && multimer == 1) {
            monoisotopicMass = (mz + adductMass) * charge;
            //System.out.println(" -M = (" + mz + " + " + adductMass + ")* " + charge + " = " + monoisotopicMass);
        }
        // CASE 3: Multimers (e.g., [2M+H]+), assume single charge
        // If Adduct is a dimer or multimer the formula is M =  (mz +- adductMass) / numberOfMultimer
        else if (charge == 1 && multimer > 1) {
            monoisotopicMass = (mz + adductMass) / multimer;
            //System.out.println(" -M = (" + mz + " + " + adductMass + ")/ " + multimer + " = " + monoisotopicMass);
        }
        // CASE 4: Multimers with multiple charges (e.g., [2M+2H]2+)
        else {
            monoisotopicMass = ((mz + adductMass) * charge) / multimer;
        }
        return monoisotopicMass;
    }

    /**
     * Calculate the mz of a monoisotopic mass with the corresponding adduct
     * Supports single and multiple charges, as well as multimers (e.g., [2M+Na]+).
     *
     * @param monoisotopicMass The neutral mass of the molecule (M).
     * @param adduct The adduct notation (e.g. "[M+H]+", "[2M+Na]+", "[M+2H]2+").
     *
     * @return The expected m/z for the adduct.
     */
    public static Double getMZFromMonoisotopicMass(Double monoisotopicMass, String adduct) {

        Double mz = null;
        Double adductMass = getAdductShift(adduct);
        // !! TODO METHOD
        // !! TODO Create the necessary regex to obtain the multimer (number before the M) and the charge (number before the + or - (if no number, the charge is 1).
        int multimer = extractMultimer(adduct);
        int charge = extractCharge(adduct);

        // Case 1: Single charge, no multimer (m/z = M +- adductMass).
        if (charge == 1 && multimer == 1) {
            mz = monoisotopicMass - adductMass;
        }
        // Case 2: Multiple charges, no multimer (mz = M/charge +- adductMass).
        else if (charge > 1 && multimer == 1) {
            mz = (monoisotopicMass / charge) - adductMass;
        }
        // Case 3: Multimer, single charge (mz = M * numberOfMultimer +- adductMass).
        else if (charge == 1 && multimer > 1) {
            mz = (monoisotopicMass * multimer) - adductMass;
        }
        // Case 4: Multimer with multiple charges (mz = (M * numberOfMultimer +- adductMass) / charge).
        else {
            mz = ((monoisotopicMass * multimer)/charge) - adductMass;
        }
        return mz;
    }

    /**
     * Returns the ppm difference between measured mass and theoretical mass
     *
     * @param experimentalMass    Mass measured by MS
     * @param theoreticalMass Theoretical mass of the compound
     */
    public static int calculatePPMIncrement(Double experimentalMass, Double theoreticalMass) {
        int ppmIncrement;
        ppmIncrement = (int) Math.round(Math.abs((experimentalMass - theoreticalMass)*1000000 / theoreticalMass));
        return ppmIncrement;
    }

    /**
     * Returns the ppm difference between measured mass and theoretical mass
     *
     * @param experimentalMass    Mass measured by MS
     * @param ppm ppm of tolerance
     */
    public static double calculateDeltaPPM(Double experimentalMass, int ppm) {
        double deltaPPM;
        deltaPPM = Math.abs((experimentalMass * ppm) / 1000000);
        return deltaPPM;
    }

    /**
     * Extracts the multimer number from the adduct string using a regex pattern.
     * The multimer refers to the number of molecules in the adduct. If no number is present,
     * the method defaults to 1.
     *
     * @param adduct The adduct string in standard format (e.g. "[2M+H]+").
     * @return The multimer count as an integer (default 1 if not found).
     */
    private static int extractMultimer(String adduct) {
        // Pattern looks for a number before 'M', e.g. [2M or [M
        Matcher m = Pattern.compile("\\[([0-9]*)M").matcher(adduct);
        if (m.find()) {
            String num = m.group(1); // Group 1 captures the digits before 'M'
            return num.isEmpty() ? 1 : Integer.parseInt(num);
        }
        return 1; // Default: no multimer specified
    }

    /**
     * Extracts the charge value from the adduct string using a regex pattern.
     * The charge is found as a number directly before the '+' or '-' at the end of the adduct.
     * If no number is present, the method assumes a default charge of 1.
     *
     * @param adduct The adduct string in standard format (e.g. "[M+2H]2+").
     * @return The charge value (default 1 if not explicitly stated).
     */
    private static int extractCharge(String adduct) {
        // Match the last digit(s) followed by + or − before the closing bracket
        Matcher m = Pattern.compile("([0-9]*)([+-])\\]?$").matcher(adduct);
        if (m.find()) {
            String num = m.group(1);  // May be empty
            return num.isEmpty() ? 1 : Integer.parseInt(num);
        }
        return 1; // Default if no explicit charge found
    }

    /**
     * Retrieves the mass shift value associated with the given adduct.
     * The method checks both the positive and negative adduct maps. If the adduct is not
     * found, it throws an exception.
     *
     * @param adduct The adduct string (e.g. "[M+H]+", "[M-H]−").
     * @return The corresponding mass shift as a double value.
     * @throws IllegalArgumentException if the adduct is not found in either map.
     */
    public static double getAdductShift(String adduct) {
        // Check in positive and negative adducts
        Map<String, Double> pos = AdductList.MAPMZPOSITIVEADDUCTS;
        Map<String, Double> neg = AdductList.MAPMZNEGATIVEADDUCTS;

        // Return the value if found
        if (pos.containsKey(adduct)) return pos.get(adduct);
        if (neg.containsKey(adduct)) return neg.get(adduct);

        // Throw error if adduct is unknown
        throw new IllegalArgumentException("Adduct not found: " + adduct);
    }
/*
    public static void main(String[] args){

        // Example for [M+H]+
        Double M1 = 700.0;
        String a1 = "[M+H]+";
        Double mz_result = getMZFromMonoisotopicMass(M1, a1);   // → ≈ 701.007
        System.out.println("Adduct " + a1 + ": " + "\n - Monoisotopic Mass: " + M1 + "\n - m/z: " + mz_result);

        Double mz1 = 701.007276;
        Double M1_result = Adduct.getMonoisotopicMassFromMZ(mz1, a1);
        System.out.println("Adduct " + a1 + ": " + "\n - m/z: " + mz1 + "\n - Monoisotopic Mass: " + M1_result);


        // Example for [M+2H]2+
        Double M2 = 700.0;
        String a2 = "[M+2H]2+";
        Double mz2_result = Adduct.getMZFromMonoisotopicMass(M2, a2);   // → ≈ 351.007
        System.out.println("Adduct " + a2 + ": " + "\n - Monoisotopic Mass: " + M2 + "\n - m/z: " + mz2_result);

        Double mz2 = 351.007276;
        Double M2_result = Adduct.getMonoisotopicMassFromMZ(mz2, a2);   // → ≈ 700.0
        System.out.println("Adduct " + a2 + ": " + "\n - m/z: " + mz2 + "\n - Monoisotopic Mass: " + M2_result);


        // Example for [2M+Na]+
        Double M3 = 700.0;
        String a3 = "[2M+Na]+";
        Double mz3_result = Adduct.getMZFromMonoisotopicMass(M3, a3);   // → ≈ 1422.989
        System.out.println("Adduct " + a3 + ": " + "\n - Monoisotopic Mass: " + M3 + "\n - m/z: " + mz3_result);

        Double mz3 = 1422.989218;
        Double M3_result = Adduct.getMonoisotopicMassFromMZ(mz3, a3);   // → ≈ 700.0
        System.out.println("Adduct " + a3 + ": " + "\n - m/z: " + mz3 + "\n - Monoisotopic Mass: " + M3_result);
    }*/

}
