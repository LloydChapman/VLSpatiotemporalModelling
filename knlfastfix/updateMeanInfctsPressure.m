function lambda_mean=updateMeanInfctsPressure(lambda_mean,lambda,i)
lambda_mean=i/(i+1)*lambda_mean+lambda/(i+1);