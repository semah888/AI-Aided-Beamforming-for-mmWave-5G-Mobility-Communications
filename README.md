# AI-Aided Beamforming for mmWave 5G Mobility Communications

**Overview**

This repository contains MATLAB and Python codes designed for AI-aided beamforming in millimeter-wave (mmWave) 5G mobility communications. The implemented applications allow for beam prediction and tracking using machine learning (ML) techniques and traditional exhaustive search approaches.

**Repository Structure**

1. `3GPP_Channel_Implementation/`: Contains MATLAB implementations for 3GPP channel models used for simulations.
2. `Traditional Beamforming/`: Implements conventional beamforming methods using DFT codebooks and exhaustive search.
3. `Computer Vision - based Beam Tracking/`: Includes Python-based object detection models (e.g., YOLOv8) for tracking vehicles to aid beam tracking.
4. `ML-driven Beam Prediction/`: Contains machine learning models (DNN, KNN, Decision Tree) for predicting optimal beams based on system parameters.
5. README.md: Documentation file describing the repository contents and functionality.

**Channel Generation and Beamforming Applications**

Two key MATLAB applications are included in the repository:
1. Channel Generator Application
- Takes network parameters as input:
  - Positions of base stations (BS) and user equipment (UE)
  - Speed of the UE
  - Antenna configuration
  - Operating frequency
  - Number of time slots
  - Update period
- Outputs:
  - Channel matrix at each time slot
  - Impulse response of the channel
  - Spectral Efficiency (SE)
  - Channel gain

2. Beamforming Application
- Takes system parameters and DFT codebook size as input
- Animates the selected beams at each time slot
- Evaluates beamforming performance using:
    - Exhaustive Search: Computes the best beam at each time slot through brute-force search.
    - ML-Driven Beamforming: Predicts beams using trained models (DNN, KNN, Decision Tree).
- Output Metrics:
    - Selected beam at the BS
    - Selected beam at the UE
    - Spectral Efficiency (SE)
    - Execution time

**Evaluation and Results**
- The repository includes MATLAB and Python scripts for evaluating the performance of traditional and AI-driven beamforming methods.
- Results demonstrate improvements in execution time and beam selection accuracy using ML models.

**Requirements**
- MATLAB (for running beamforming and channel modeling applications)
- Python (for computer vision-based beam tracking and ML-based beam prediction)
- Required Python libraries: numpy, pandas, seaborn, matplotlib, scikit-learn, keras, cv2 (for YOLOv8)

**Usage**
- Channel Generation: Run the MATLAB script in 3GPP_Channel_Implementation/ to generate the channel response.
- Beamforming: Use the MATLAB scripts in Traditional Beamforming/ or ML-driven Beam Prediction/ to perform beam selection.
- Evaluation: Analyze results using provided scripts for SE and execution time comparisons.

**Acknowledgments**
- This project is based on 3GPP standards and leverages machine learning for optimizing beamforming in 5G mobility scenarios.
