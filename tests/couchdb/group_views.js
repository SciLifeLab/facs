//Group by Ecoli bloom
function(doc) {
    if(doc.bloom_filter && doc.sample) {
        var bloom_file = doc.bloom_filter.replace(/^.*[\\\/]/, '');
        var sample_file = doc.sample.replace(/^.*[\\\/]/, '');
        if(bloom_file == 'eschColi_k12.bloom') {
            emit([bloom_file, sample_file], doc._id);
        }
    }
}

//Group by Dm3 bloom
function(doc) {
    if(doc.bloom_filter && doc.sample) {
        var bloom_file = doc.bloom_filter.replace(/^.*[\\\/]/, '');
        var sample_file = doc.sample.replace(/^.*[\\\/]/, '');
        if(bloom_file == 'dm3.bloom') {
            emit([bloom_file, sample_file], doc._id);
        }
    }
}

//Group by human bloom
function(doc) {
    if(doc.bloom_filter && doc.sample) {
        var bloom_file = doc.bloom_filter.replace(/^.*[\\\/]/, '');
        var sample_file = doc.sample.replace(/^.*[\\\/]/, '');
        if(bloom_file == 'hg19.bloom') {
            emit([bloom_file, sample_file], doc._id);
        }
    }
}

//Group by phiX bloom
function(doc) {
    if(doc.bloom_filter && doc.sample) {
        var bloom_file = doc.bloom_filter.replace(/^.*[\\\/]/, '');
        var sample_file = doc.sample.replace(/^.*[\\\/]/, '');
        if(bloom_file == 'phiX.bloom') {
            emit([bloom_file, sample_file], doc._id);
        }
    }
}

