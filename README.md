Hyperlink Prediction Toolbox
============================

About
-----

A hyperlink relaxes the restriction in traditional link prediction that two nodes form a link. Instead, it allows an arbitrary number of nodes to jointly form a multiway relation. Hyperlink prediction is much harder than link prediction, as the total number of hyperlinks is O(2^n) instead of O(n^2). 

This toolbox contains the source code of our paper:

> M. Zhang,  Z. Cui,  S. Jiang,  and Y. Chen,  Beyond Link Prediction: Predicting Hyperlinks in Adjacency Space,  Proc. AAAI Conference on Artificial Intelligence (AAAI-18). [\[PDF\]](http://www.cse.wustl.edu/~muhan/papers/AAAI_2018_Hyperlink.pdf)

which studies the hyperlink prediction problem. We proposed a novel Coordinated Matrix Minimization (CMM) algorithm, which predicts hyperlinks in the vertex adjacency space and achieves state-of-the-art performance on two tasks: recipe prediction and metabolic reaction prediction. This toolbox also includes implementations of various hyperlink prediction baselines.

Have a try by typing "Main_meta" in MATLAB!

How to run
----------

The file "Main_meta.m" is the main program for running the metabolic reaction prediction experiments. The file "Main_dish" is for dish prediction. Please change experimental settings inside these two files.

The folder "data/" contains the preprocessed data used in the experiments. They are: iAB_RBC_283.mat, iAF1260b.mat, iAF692.mat, iHN637.mat, iIT341.mat and iJO1366.mat (each is a metabolic network model), and chuancai.txt, yuecai.txt (collections of Sichuan and Cantonese recipes). The folder "utils/" contains some preprocessing scripts, which may be useful if you want to experiment on more datasets.

By default, "Main_meta.m" will save results in "result/". There is a "Meta_Plot.m" in "utils/" which plots figures according to the results.

Required libraries
------------------

The software symnmf-master is needed and already included in "software/". 

The software liblinear is required for some baselines. We have included a "liblinear.tar.gz" in "software/". Type:

    tar -xvzf liblinear.tar.gz

to unzip it. Then please follow its README to install it (and its MATLAB version, see README in /matlab/). After installation, you need to change the names of "train.mexa64" and "predict.mexa64" to "liblinear_train.mexa64" and "liblinear_predict.mexa64" respectively, otherwise they will cause conflictions.

The software libFM is required for some baselines. We have included a "libfm-1.42.src.tar.gz" in "software/". Type:

    tar -xvzf libfm-1.42.src.tar.gz

to unzip it. Then, compile it by

    cd libfm-1.42-src
    make all
    cp bin/libFM ../../../data/FM_temp/

This will compile and copy the executable libFM into "data/FM_temp".

Reference
---------

If you find the code useful, please cite our paper:

    @paper{AAAI1817136,
        author = {Muhan Zhang and Zhicheng Cui and Shali Jiang and Yixin Chen},
        title = {Beyond Link Prediction: Predicting Hyperlinks in Adjacency Space},
        conference = {AAAI Conference on Artificial Intelligence},
        year = {2018},
        keywords = {link prediction; hypergraph; hyperlink; EM; metabolic network},
        url = {https://www.aaai.org/ocs/index.php/AAAI/AAAI18/paper/view/17136/16754}
    }

Muhan Zhang, muhan@wustl.edu
12/2/2017
