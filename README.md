# Speech-Digit-Recognition
Feature Extraction

    Signal Preprocessing:
        DC Offset Correction: Removes any constant offset from the signal.
        Normalization: Adjusts the amplitude of the signal for uniformity.
        Frame Selection: Extracts steady-state frames for further analysis.

    Frame-Level Processing:
        Hamming Window Application: Smooths the edges of each frame.
        Autocorrelation: Computes a measure of similarity of the signal with a delayed version of itself.
        LPC and Cepstral Coefficients:
            LPC (Linear Predictive Coding) coefficients are calculated.
            Cepstral coefficients are derived for each frame and weighted for emphasis.

    Output Generation:
        Writes the weighted cepstral coefficients for each frame into a structured file. and after that merged all file and make a universe.csv file




Clustering (LBG Algorithm)

    Initialization:
        Loads a "universe" of feature vectors from UNIVERSE.csv.
        Begins with one centroid (average of all vectors).

    Iterative Clustering:
        Splits centroids iteratively, doubling the codebook size.
        Uses the Tokhura distance to assign vectors to clusters and update centroids.
        Stops when the distortion stabilizes within a threshold.

    Codebook Creation:
        Outputs a codebook (codebook.csv) with 32 vectors, each having 12 dimensions.



Observation Sequence Generation

    Mapping to Nearest Codebook Vector:
        Reads the feature vector file for each sample.
        For each frame, finds the closest codebook vector using the Tokhura distance.
        Outputs the sequence of observations (OUTPUT_OBSERVATION files).


The training phase (step 4) initializes and calculates parameters such as the transition matrix A B for every digit separately and it is stored in different file as well and finally average of all value is taken to generate the final model.


Then testing is done using HMM1 as on which out of the 10 model it is more likely to be generated and output is given.

Then searching of sentence or text in the file paragraph.txt is written 





If you want to run for different digit (like i have run for code or digit 1 in the video ) just change the input parameter in the answer function which is present in the startrecord function(answer function is present in last line of it and start reord function is present in just before 2 function below main function) to the digit you are going to speak. It is basically to check wheather it is recognising the code correctly or not. We can also remove it but I have kept it check if it is recognizing correctly or not.
