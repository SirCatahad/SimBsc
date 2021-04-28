# SimBsc
The description in the original paper seems weird, the r-squared values of 0, .99 and .5 do not match the (respective) covariance values of 0, 3.33 and 10 as described in the paper
3.33 would correspond to .1, better to .11
.5 does not correspond to 10 (.99 does)

Since later in the paper they refer to "0, .1 and .5 again" I assume that .1 is meant for covariance of 3.33 (although it's more like .11) and perhaps .99 was meant for covariance of 10, although I might be completely off.
   
I have decided not to use r.squared values but instead keep the covariances fixed for the second study due to the issues mentioned above. If it is a misunderstanding it would be a simple fix to replace covariance by r.squared in the data generation process
