function error = getHrCurlNormforProblem3(c,globalA)

    error = c' * globalA * c;
    error = sqrt(error);

end