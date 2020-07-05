# Reinforcement Ductile Fracture Simulation

Recent research has used empirical models (e.g., Manson, 1954; Coffin, 1954) to evaluate the cyclic strain lives of reinforcing steel materials (e.g., Mander et al., 1994; Brown and Kunnath, 2000; Aoyama, 2001; Slavin and Ghannoum, 2015; Zhong and Deierlein, 2019).  Typically, for a specific type of reinforcement, constant-amplitude low-cycle fatigue tests were conducted to measure half-cycle numbers under different strain amplitudes. Then modeling coefficients were calibrated, based on the assumption that the half-cycle number and the strain amplitude has a log-log linear relation.  

However, one common disadvantage of empirical models is that the prediction accuracy is not always guaranteed by the log-log linearity under large strain amplitudes (e.g., strain and stress concentrations in the necking stage).  Researchers also agree that the actual failure process after bar buckling is more complex since localized residual strains accelerate the material failure (e.g., Sokoli, 2018; Barcley and Kowalsky, 2019).  This buckling effect, which is also dependent on the loading path of the reinforcement, is not captured by empirical damage models.  Consequently, in a strict sense, these empirical models can help estimate the probability of the first bar fracture in concrete structures but are incapable of explicitly simulating reinforcement fractures.

<p align="center">
 <img width="80%" height="80%" src="https://github.com/kuanshi/ductile-fracture/blob/master/doc/image/RDFM_issue.png">
</p>

The Reinforcement Ductile Fracture Model (RDFM) is developed (1) to implement the physical mechanism of void-growth (Kanvinde and Deierlein, 2006, Kanvinde and Deierlein, 2007, Fell et al., 2009, Myers et al., 2010, Smith et al., 2017; and Masao, 2018) for more accurate predictions under monotonic and ultra-low cyclic loads, (2) to explicitly track local strain concentration after bar buckling, (3) to provide a [probabilistic description of bar fracture](https://github.com/kuanshi/ductile-fracture/blob/master/doc/image/RDFM_issue.png), and (4) to enable [explicit bar fracture simulations](https://github.com/kuanshi/ductile-fracture/blob/master/doc/image/RDFM_fracture.png) in [reinforced concrete beam-column fiber elements](https://github.com/kuanshi/ductile-fracture/blob/master/doc/image/RDFM_model.png). The figure below illustrates the workflow of the RDFM.

<p align="center">
 <img width="60%" height="60%" src="https://github.com/kuanshi/ductile-fracture/blob/master/doc/image/RDFM.png">
</p>
