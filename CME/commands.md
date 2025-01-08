# Terminal Commands to Run Jobs on Delta


## Login into Delta login node

```bash
ssh USERNAME@login.delta.ncsa.illinois.edu
```

***Replace*** `USERNAME` with your Delta username. You need to type your password and do 2FA.

##  Enter your projects directory and copy CME tutorials

Navigate to your directory

```bash
cd /projects/bddt/$USER
```

Copy the prepared materials to your directory. This step may take several minutes since the large Apptainer file and prepared CME/ODE Whole-Cell Model trajactories.

```bash
cp -r /projects/bddt/LM ./
```

## Launch Jupyter Notebook on Delta
You will use Jupyter Notebook to run Tutorial 1, 2 and the analysis part of Tutorial 3. The advantage of Jupyter Notebook is that you could navigate the folders and run the *.ipynb* file using GUI.  
+ **First**: Submit a job to a Delta GPU node.  
     Here *srun* launch interactive job onto Delta, *partition* claims A100 GPU, and for 4 hours *time*. We need to specify the `Port` for *jupyter-notebook*.  
    ***Copy*** the following command and ***Replace*** `Port` with a four digit non-trivial number to avoid using the same port as others. Group A should use 1111, 2222, 3333 or 4444. Group B can use 5555, 6666, 7777 or 9999.
    ```bash
    srun --account=bddt-delta-gpu --partition=gpuA100x4 --time=04:00:00 --mem=64g --gpus-per-node=1 --tasks-per-node=1 --cpus-per-task=16 --nodes=1 apptainer exec --nv --containall --bind /projects/bddt/$USER/:/workspace /projects/bddt/$USER/LM/LM.sif jupyter-notebook /workspace/ --no-browser --port=Port --ip=0.0.0.0 --allow-root
    ```  
    Then you should wait for Delta to allocate the resources for you that usually takes less than 1 minute. When you see similar things as the following, you are good to proceed:
    ```bash
    srun: job 3546627 queued and waiting for resources
    srun: job 3546627 has been allocated resources
    WARNING: could not mount /etc/localtime: not a directory
    [I 19:07:57.203 NotebookApp] Writing notebook server cookie secret to /u/$USER/.local/share/jupyter/runtime/notebook_cookie_secret
    [I 19:07:58.314 NotebookApp] [jupyter_nbextensions_configurator] enabled 0.6.3
    [I 19:07:58.316 NotebookApp] Serving notebooks from local directory: /workspace
    [I 19:07:58.316 NotebookApp] Jupyter Notebook 6.4.12 is running at:
    [I 19:07:58.316 NotebookApp] http://`DeltaNode`.ncsa.illinois.edu:8811/?token=b2e7ca15cd9dc3a6893a1273e359c88869225bc29d66c80c
    [I 19:07:58.316 NotebookApp]  or http://127.0.0.1:$Port/?token=b2e7ca15cd9dc3a6893a1273e359c88869225bc29d66c80c
    [I 19:07:58.316 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
    [C 19:07:58.329 NotebookApp]

        To access the notebook, open this file in a browser:
            file:///u/$USERNAME/.local/share/jupyter/runtime/nbserver-13-open.html
        Or copy and paste one of these URLs:
            http://`DeltaNode`.delta.ncsa.illinois.edu:$Port/?token=b2e7ca15cd9dc3a6893a1273e359c88869225bc29d66c80c
        or http://127.0.0.1:$Port/?token=b2e7ca15cd9dc3a6893a1273e359c88869225bc29d66c80c
    ```

    The last two line contains the Delta GPU node `DeltaNode`, which is the node where your job is running on.

+ **Second**: SSH into the Delta GPU node.  
  Open a second terminal and Copy the following command.  
  Your `DeltaNode` can be found from the information above in last two lines after `http://`. ***Replace*** `DeltaNode` with your node you see above and ***Replace*** `USERNAME` with your Delta username. ***Replace*** two `Port` with the 4 digit number you used.
    
    ```bash
    ssh -l USERNAME  -L 127.0.0.1:Port:DeltaNode.delta.internal.ncsa.edu:Port dt-login.delta.ncsa.illinois.edu
    ```

    You need to type you password and do 2FA AGAIN.

+ **Third**: Open Jupyter Notebook Webpage.   
  Copy the last URL in the first terminal and paste to one browser (Firefox, Chrome, ...) to open Jupyter Notebook.


## Run CME/ODE Whole-Cell Model's Python Scripts

Run this part when you start Tutorial 3, CME/ODE Whole-Cell Model of Minimal Cell.

You will submit a job to run Whole-Cell Model in parallel on Delta GPU node.

### Open a new terminal and then login to Delta ###

+ **First**: Go to *programs* folder

    ``` bash
    cd /projects/bddt/$USER/LM/CME/WholeCellModel/programs
    ```

+ **Second**: Submit the bash file

    ```bash
    sbatch mpirun.sh
    ```
     In the given bash file, you will launch 2 minutes simulation of 4 replicates. You are encouraged to change the parameters to run longer of more replicates. Normally, it takes 6 physical hours to finish 2 biological hours' simulation on a Delta A100 GPU node.
  
+ **Third**: Check the status of your job.  
    *PD* means waiting to run, *R* running.

    ```bash
    squeue -u $USER
    ```
    go to *output_4replicates* folder 
    
    ``` bash
    cd /projects/bddt/$USER/LM/CME/WholeCellModel/output_4replicates
    ```

### Run pkl.ipynb then plotting.ipynb on Jupyter Notebook Webpage
+ **First**: Navigate to `./LM/CME/WholeCellModel/analysis` in Jupyter Notebook Webpage.
  
+ **Second**: Open *pkl.ipynb* and **Run All** to serialize the CSV files into one single python pickle file.  
  10 CSV files are generated beforehand for the analysis.
  
+ **Thrid**: Open *plotting.ipynb* and **Run All** to generate the plots.  
  Compare the plots with figures in *Summer_School_CME_2024.pdf*.

