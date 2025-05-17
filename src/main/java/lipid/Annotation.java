package lipid;

import adduct.Adduct;
import adduct.AdductList;

import java.util.*;

/**
 * Class to represent the annotation over a lipid
 */
public class Annotation {

    private final Lipid lipid;
    private final double mz;
    private final double intensity; // intensity of the most abundant peak in the groupedPeaks
    private final double rtMin;
    private final IoniationMode ionizationMode;
    private String adduct; // !!TODO The adduct will be detected based on the groupedSignals
    private final Set<Peak> groupedSignals;
    private int score;
    private int totalScoresApplied;


    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param ionizationMode
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IoniationMode ionizationMode) {
        this(lipid, mz, intensity, retentionTime, ionizationMode, Collections.emptySet());
    }

    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param ionizationMode
     * @param groupedSignals
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IoniationMode ionizationMode, Set<Peak> groupedSignals) {
        this.lipid = lipid;
        this.mz = mz;
        this.rtMin = retentionTime;
        this.intensity = intensity;
        this.ionizationMode = ionizationMode;
        // !!TODO This set should be sorted according to help the program to deisotope the signals plus detect the adduct
        this.groupedSignals = new TreeSet<>(Comparator.comparing(Peak::getMz));
        this.groupedSignals.addAll(groupedSignals);
        this.score = 0;
        this.totalScoresApplied = 0;
        detectAdductFromPeaks();
    }

    public Lipid getLipid() {
        return lipid;
    }

    public double getMz() {
        return mz;
    }

    public double getRtMin() {
        return rtMin;
    }

    public String getAdduct() {
        return adduct;
    }

    public void setAdduct(String adduct) {
        this.adduct = adduct;
    }

    public double getIntensity() {
        return intensity;
    }

    public IoniationMode getIonizationMode() {
        return ionizationMode;
    }

    public Set<Peak> getGroupedSignals() {
        return Collections.unmodifiableSet(groupedSignals);
    }


    public int getScore() {
        return score;
    }

    public void setScore(int score) {
        this.score = score;
    }

    // !CHECK Take into account that the score should be normalized between -1 and 1
    public void addScore(int delta) {
        this.score += delta;
        this.totalScoresApplied++;
    }

    /**
     * @return The normalized score between 0 and 1 that consists on the final number divided into the times that the rule
     * has been applied.
     */
    public double getNormalizedScore() {
        return (double) this.score / this.totalScoresApplied;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Annotation)) return false;
        Annotation that = (Annotation) o;
        return Double.compare(that.mz, mz) == 0 &&
                Double.compare(that.rtMin, rtMin) == 0 &&
                Objects.equals(lipid, that.lipid);
    }

    @Override
    public int hashCode() {
        return Objects.hash(lipid, mz, rtMin);
    }

    @Override
    public String toString() {
        return String.format("Annotation(%s, mz=%.4f, RT=%.2f, adduct=%s, intensity=%.1f, score=%d)",
                lipid.getName(), mz, rtMin, adduct, intensity, score);
    }

    /**
     * Asignar un aducto a la variable local a partir de los peaks de la lista de adducts
     * usando la mz de la anotación como referencia para determinar el aducto
     */
    public void detectAdductFromPeaks() {
        List<Peak> peaks = new ArrayList<>(groupedSignals);
        double annotationMz = this.mz; // mz de la anotación que usar como referencia

        // Primero intentamos con los aductos positivos
        for (Map.Entry<String, Double> adduct : AdductList.MAPMZPOSITIVEADDUCTS.entrySet()) {
            // Calculamos la masa monoisotópica de la anotación con este aducto
            double monoisotopicMass = Adduct.getMonoisotopicMassFromMZ(annotationMz, adduct);

            // Si encontramos picos que correspondan a esta masa, asignamos este aducto
            if (matchPeaksWithMonoisotopicMass(peaks, monoisotopicMass, adduct.getKey())) {
                return;
            }
        }

        // Si no encontramos con los positivos, intentamos con los negativos
        for (Map.Entry<String, Double> adduct : AdductList.MAPMZNEGATIVEADDUCTS.entrySet()) {
            // Calculamos la masa monoisotópica de la anotación con este aducto
            double monoisotopicMass = Adduct.getMonoisotopicMassFromMZ(annotationMz, adduct);

            // Si encontramos picos que correspondan a esta masa, asignamos este aducto
            if (matchPeaksWithMonoisotopicMass(peaks, monoisotopicMass, adduct.getKey())) {
                return;
            }
        }
    }

    /**
     * Verifica si alguno de los picos corresponde a la masa monoisotópica calculada
     * @param peaks lista de picos a verificar
     * @param monoisotopicMass masa monoisotópica calculada
     * @param adductKey clave del aducto
     * @return true si se encontró correspondencia
     */
    private boolean matchPeaksWithMonoisotopicMass(List<Peak> peaks, double monoisotopicMass, String adductKey) {
        int PPM_TOLERANCE = 10;
        double deltaPPM = Adduct.calculateDeltaPPM(monoisotopicMass, PPM_TOLERANCE);

        for (Peak peak : peaks) {
            // Para cada pico, probamos todos los aductos posibles
            for (Map.Entry<String, Double> peakAdduct : AdductList.MAPMZPOSITIVEADDUCTS.entrySet()) {
                double peakMonoisotopicMass = Adduct.getMonoisotopicMassFromMZ(peak.getMz(), peakAdduct);

                if (Math.abs(monoisotopicMass - peakMonoisotopicMass) <= deltaPPM) {
                    // Asignamos el aducto de la anotación
                    this.adduct = adductKey;

                    System.out.println("Relación de aductos encontrada:");
                    System.out.println("  Anotación: " + this.mz + " con aducto " + adductKey);
                    System.out.println("  Peak: " + peak.getMz() + " con aducto " + peakAdduct.getKey());
                    System.out.println("  Masas monoisotópicas: " + monoisotopicMass + " ~ " + peakMonoisotopicMass +
                            " (ΔDa: " + Math.abs(monoisotopicMass - peakMonoisotopicMass) + " ≤ " + deltaPPM + ")");

                    return true;
                }
            }

            // También probamos con los aductos negativos
            for (Map.Entry<String, Double> peakAdduct : AdductList.MAPMZNEGATIVEADDUCTS.entrySet()) {
                double peakMonoisotopicMass = Adduct.getMonoisotopicMassFromMZ(peak.getMz(), peakAdduct);

                if (Math.abs(monoisotopicMass - peakMonoisotopicMass) <= deltaPPM) {
                    // Asignamos el aducto de la anotación
                    this.adduct = adductKey;

                    System.out.println("Relación de aductos encontrada:");
                    System.out.println("  Anotación: " + this.mz + " con aducto " + adductKey);
                    System.out.println("  Peak: " + peak.getMz() + " con aducto " + peakAdduct.getKey());
                    System.out.println("  Masas monoisotópicas: " + monoisotopicMass + " ~ " + peakMonoisotopicMass +
                            " (ΔDa: " + Math.abs(monoisotopicMass - peakMonoisotopicMass) + " ≤ " + deltaPPM + ")");

                    return true;
                }
            }
        }

        return false;
    }

}
