function batch = assign_batch_id(batch, model)
% ASSIGN_BATCH_ID - Assigns a unique identifier to a batch based on a model name and random hash.
%
% This function generates a unique identifier for a batch by combining the model name with a randomly generated hash code. 
% The identifier is stored in the batch struct under the 'id' field, enabling tracking or referencing of the batch in subsequent 
% operations, such as data storage or retrieval.
%
% INPUTS:
%   batch      - A struct representing the batch, which will be updated with an 'id' field.
%   model      - A struct containing at least a 'name' field, representing the model name (string).
%
% OUTPUTS:
%   batch      - The updated batch struct with a new 'id' field:
%                - batch.id: A string combining the model name and a random hash code (e.g., "ModelName1234567890").
%
% FUNCTIONALITY:
% - Generates a random integer hash code between 1 billion (1e9) and 10 billion minus 1 (1e10-1).
% - Concatenates the model name with the hash code to create a unique identifier.
% - Assigns the resulting string to the 'id' field of the batch struct.
%
% NOTES:
% - The hash code is a large random integer to minimize the likelihood of duplicate IDs across batches.
% - The model name is assumed to be a string stored in `model.name`.
% - No validation is performed to ensure uniqueness; the large range of possible hash codes makes collisions unlikely.
%
% Author: Ludwig Horvath
% Date: 2/11/2025


hashcode = randi([1e9, 1e10-1]);

name = append(model.name, string(hashcode));

batch.id = name;

end