function(doc) {
  if (doc.sample && doc.contamination_rate) {
    // only check organisms screened against ecoli for now
    var bloom_file = doc.bloom_filter.replace(/^.*[\\\/]/, '')
    var sample_file = doc.sample.replace(/^.*[\\\/]/, '')
    if (bloom_file === 'eschColi_K12.bloom') {
       emit([sample_file, bloom_file], doc.contamination_rate);
    }
  }
}
