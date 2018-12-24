public class sgpfunction {
    public static void go(String satelliteNumber, String internationalDesignator, int epochYear, int epoch, int firstTimeDerivativeofMeanMotion, int secondTimeDerivativeofMeanMotion, int bstar, int ephemerisType, int elementNumber, int checksum, int inclination, int rightAscensoinoftheAscendingNode, int eccentricity, int argumentOfPerigee, int meanAnomaly, int meanMotion, int revolutionNumberAtEpoch) {
        System.out.println("1 " + satelliteNumber + "U " + internationalDesignator + " " + Integer.toString(epochYear) + Integer.toString(epoch) + Integer.toString(firstTimeDerivativeofMeanMotion) + " " + Integer.toString(bstar) + " " + Integer.toString(ephemerisType) + " " + Integer.toString(elementNumber) + " " +  Integer.toString(checksum));
        System.out.println("2 " + satelliteNumber + " " + Integer.toString(inclination) + " " + Integer.toString(rightAscensoinoftheAscendingNode) + " " + Integer.toString(eccentricity) + " " + Integer.toString(argumentOfPerigee) + " " + Integer.toString(meanAnomaly) + " " + Integer.toString(meanMotion) + " " + Integer.toString(revolutionNumberAtEpoch) + Integer.toString(checksum));
    }
}
