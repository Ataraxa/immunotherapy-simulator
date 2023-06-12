function error = quick_testing(input_struct)
    error = abs(input_struct.v1 - input_struct.v2 + 10);
end