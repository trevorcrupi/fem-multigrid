function error = getHrCurlErrorforProblem3(c,globalA)

    error = c' * globalA * c;
    error = sqrt(error);

end