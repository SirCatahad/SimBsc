# SimBsc
Current issues:
 - Data generation seems to work fine
 - PDMI (mice "norm") seems wonky from the results. Not sure what exactly is wrong
 - MCCI (Monte Carlo Confidence Intervals) not correctly implemented (pretty sure about that)
 - Might also be Monte Carlo Standard Errors.
 - Coverage of CIs seems way off, 50-60% of CIs cover true parameter seems unlikely in 100+ iterations
 - GGPlot does not plot when running everything. Workaround is to rerun any command so it can throw it's timeout, then manually rerun the plot commands.
   Might be a bug in the package itself.
   
