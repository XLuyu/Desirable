package edu.nus

import java.io.File

class KmerAnalyzer(val settings: Desirable) {
    val bf = BloomFilter(settings.K,8*3000000000L+7, listOf(2,3,5,7,11,13))
    val troughAndPeak = settings.countFiles.map { getSpectrumTroughPeak(it) }

    private fun generateKmerSpectrum(countFile: File): HashMap<Int, Long> {
        val spectrum = HashMap<Int, Long>()
        countFile.forEachLine {
            val freq = it.split("\t")[1].toInt()
            spectrum[freq] = spectrum.getOrDefault(freq, 0) + 1
        }
        return spectrum
    }
    fun getSpectrumTroughPeak(countFile: File): Pair<Int, Int> { // get first local minimum in kmer spectrum
        val spectrum = generateKmerSpectrum(countFile).toList().sortedBy { it.first }
        val total = spectrum.map { it.first * it.second }.sum()
        var accumulate = 0L
        var trough = spectrum.indices.first {
            if (spectrum[it].second <= spectrum[it + 1].second) return@first true
            accumulate += spectrum[it].first * spectrum[it].second
            (accumulate > total * 0.5)
        }
        if (accumulate > total * 0.5) trough = 0
        val peak = spectrum.filter { it.first>=trough }.maxBy { it.first * it.second }!!.first
        println("[KmerSpectrum]\nTrough (Error threshold): $trough\nPeak (Coverage): $peak")
        return Pair(trough, peak)
    }

    fun filterTargetKmerByControl(): HashMap<Long,Long> {
        val filesReader = settings.countFiles.map { it.bufferedReader() }
        val record = Array<List<String>>(filesReader.size,{ listOf(" "," ") })
        val kmers = HashMap<Long,Long>()
        while (filesReader[0].ready()) {
            record[0] = filesReader[0].readLine().split('\t')
            if (record[0][1].toInt() < troughAndPeak[0].first) continue
            bf.insert(Util.canonical(record[0][0]))
            for ( i in 1 until filesReader.size)
                while ((record[0][0] > record[i][0]) && filesReader[i].ready())
                    record[i] = filesReader[i].readLine().split('\t')
            if ((1 until record.size).any { i -> record[i][0] == record[0][0] && record[i][1].toInt()>=troughAndPeak[i].first}) continue
            kmers[Util.canonical(record[0][0])] = record[0][1].toLong()
        }
        println("[Kmer] ${kmers.size} exclusive Kmers")
        return kmers
    }
    fun filterTargetByControl(): HashMap<Long,Long>  {
        val exKmer = filterTargetKmerByControl()
        for (f in settings.files[0]){
            val seqReader = ReadFileReader(f,"Exclusive Read Screen")
            while (seqReader.nextRead()){
                val kmerList = seqReader.decomposeKmer(settings.K)
                if (!kmerList.any{ it in exKmer}) continue
                kmerList.filter { bf.contains(it) } .forEach { exKmer[it] = 1 }
            }
        }
        val targetReader = settings.countFiles[0].bufferedReader()
        while (targetReader.ready()) {
            val (targetKmer, targetFreq) = targetReader.readLine().split('\t')
            val kmer = Util.canonical(targetKmer)
            val freq = targetFreq.toLong()
            if (!exKmer.contains(kmer)) continue
            if (freq < troughAndPeak[0].first) exKmer.remove(kmer) else exKmer[kmer] = freq
        }
        return exKmer
    }
}