package edu.nus

import java.io.File

class KmerAnalyzer(val targetCount: File, val controlCount: File) {

    fun generateKmerSpectrum(countFile: File): HashMap<Int, Long> {
        val spectrum = HashMap<Int, Long>()
        countFile.forEachLine {
            val freq = it.split("\t")[1].toInt()
            spectrum[freq] = spectrum.getOrDefault(freq, 0) + 1
        }
        return spectrum
    }

    fun getSpectrumTroughPeak(countFile: File): Pair<Int, Int> { // get first local minimum in kmer spectrum
        val spectrum = generateKmerSpectrum(countFile).toList()
        val total = spectrum.map { it.first * it.second }.sum()
        var accumulate = 0L
        var trough = spectrum.indices.first {
            if (spectrum[it].second <= spectrum[it + 1].second) return@first true
            accumulate += spectrum[it].first * spectrum[it].second
            (accumulate > total * 0.5)
        }
        if (accumulate > total * 0.5) trough = 0
        val peak = spectrum.maxBy { it.first * it.second }!!.first
        println("[KmerSpectrum]\nTrough (Error threshold): $trough\nPeak (Coverage): $peak")
        return Pair(trough, peak)
    }

    fun filterTargetByControl(targetTroughPeak: Pair<Int, Int>, controlTroughPeak: Pair<Int, Int>): HashSet<String> {
        val targetReader = targetCount.bufferedReader()
        val controlReader = controlCount.bufferedReader()
        var (controlKmer, controlFreq) = Pair(" "," ")
        val kmers = HashSet<String>()
        while (targetReader.ready()) {
            val (targetKmer, targetFreq) = targetReader.readLine().split('\t')
            if (targetFreq.toInt() < targetTroughPeak.first) continue
            while ((targetKmer > controlKmer) && controlReader.ready()){
                val tokens = controlReader.readLine().split('\t')
                controlKmer = tokens[0]
                controlFreq = tokens[1]
            }
            if (targetKmer == controlKmer && controlFreq.toInt()>=controlTroughPeak.first) continue
            kmers += targetKmer
        }
        println("[FilterKmer] ${kmers.size} Kmers left")
        return kmers
    }
}