function(doc) {
  if (doc.sample && doc.contamination_rate) {
    // only check organisms screened against ecoli for now
    var bloom_file = doc.bloom_filter.replace(/^.*[\\\/]/, '');
    var sample_file = doc.sample.replace(/^.*[\\\/]/, '');
    if (bloom_file === 'eschColi_K12.bloom'){
       if (sample_file.substring(0,6) === 'simngs') {
           if (doc.timestamp) {
               // Javascript Date.parse() does not support subsecond error format
               var timestamp = doc.timestamp.replace(/\+\d{4}/, '');
               var d1 = Date.parse(timestamp);
               if(sample_file.match(/simngs_\w+_\d+\.fastq/g)){
                            emit([bloom_file, sample_file], doc._id);
               }

               //var diff = end_run - begin_run;
               emit([sample_file, bloom_file], [d1, doc.contamination_rate]);
           }
       }
    }
  }
}
