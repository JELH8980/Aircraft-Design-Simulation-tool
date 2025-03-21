function save_batch(batch)

info = batch.info;

if exist('Database', 'dir')
    if exist('Batches', 'dir')

    else
        cd('Database\')
        mkdir('Batches')
        cd('Batches')
    end
else
    mkdir('Database');
    cd('Database\')
    mkdir('Batches')
    cd('Batches')
end

save(append(string(info.id), ".mat"), "info", "-mat")

cd('..')
cd('..')

end