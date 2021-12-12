To compile code:

nvcc .\fluidSim_seq_no_GUI.cpp -o .\fluidSim_seq_no_GUI_64
nvcc .\fluidSim_seq_no_GUI.cpp -o .\fluidSim_seq_no_GUI_128
nvcc .\fluidSim_seq_no_GUI.cpp -o .\fluidSim_seq_no_GUI_256
nvcc .\fluidSim_seq_no_GUI.cpp -o .\fluidSim_seq_no_GUI_512
nvcc .\fluidSim_seq_no_GUI.cpp -o .\fluidSim_seq_no_GUI_1024


To run tests:

.\fluidSim_seq_no_GUI_128.exe 10 0;
.\fluidSim_seq_no_GUI_128.exe 20 0;
.\fluidSim_seq_no_GUI_128.exe 30 0;
.\fluidSim_seq_no_GUI_128.exe 40 0;
.\fluidSim_seq_no_GUI_128.exe 50 0;

.\fluidSim_seq_no_GUI_64.exe 10 0;
.\fluidSim_seq_no_GUI_128.exe 10 0;
.\fluidSim_seq_no_GUI_256.exe 10 0;
.\fluidSim_seq_no_GUI_512.exe 10 0;
.\fluidSim_seq_no_GUI_1024.exe 10 0;

