package edu.nus

import java.io.File

class KmerAnalyzer(val settings: Desirable) {
    val bf = BloomFilter(settings.K,8*3000000000L+7, listOf(2,3,5,7,11,13))
    fun generateKmerSpectrum(countFile: File): HashMap<Int, Long> {
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

    fun filterTargetKmerByControl(targetTroughPeak: Pair<Int, Int>, controlTroughPeak: Pair<Int, Int>): HashMap<Long,Long> {
        val targetReader = settings.targetCount!!.bufferedReader()
        val controlReader = settings.controlCount!!.bufferedReader()
        var (controlKmer, controlFreq) = Pair(" "," ")
        val kmers = HashMap<Long,Long>()
        while (targetReader.ready()) {
            val (targetKmer, targetFreq) = targetReader.readLine().split('\t')
            if (targetFreq.toInt() < targetTroughPeak.first) continue
            bf.insert(Util.canonical(targetKmer))
            while ((targetKmer > controlKmer) && controlReader.ready()){
                val tokens = controlReader.readLine().split('\t')
                controlKmer = tokens[0]
                controlFreq = tokens[1]
            }
            if (targetKmer == controlKmer && controlFreq.toInt()>=controlTroughPeak.first) continue
            kmers[Util.canonical(targetKmer)] = targetFreq.toLong()
        }
        println("[Kmer] ${kmers.size} exclusive Kmers")
        return kmers
    }
    fun filterTargetByControl(target: List<File>, targetTroughPeak: Pair<Int, Int>, controlTroughPeak: Pair<Int, Int>): HashMap<Long,Long>  {
        val exKmer = filterTargetKmerByControl(targetTroughPeak, controlTroughPeak)
        for (f in target){
            val seqReader = ReadFileReader(f,"Exclusive Read Screen")
            while (seqReader.nextRead()){
                val kmerList = seqReader.decomposeKmer(settings.K)
                if (!kmerList.any{ it in exKmer}) continue
                kmerList.filter { bf.contains(it) } .forEach { exKmer[it] = 1 }
            }
        }
        val targetReader = settings.targetCount!!.bufferedReader()
        while (targetReader.ready()) {
            val (targetKmer, targetFreq) = targetReader.readLine().split('\t')
            val kmer = Util.canonical(targetKmer)
            val freq = targetFreq.toLong()
            if (!exKmer.contains(kmer)) continue
            if (freq < targetTroughPeak.first) exKmer.remove(kmer) else exKmer[kmer] = freq
        }
        return exKmer
    }
}