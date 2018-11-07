package edu.nus

import java.io.File

class KmerOperator(val readFiles: Array<ArrayList<File>>) {
    var sampleNames = readFiles.map { it.first().absolutePath }.toTypedArray()
    var countFiles = Array<File?>(sampleNames.size,{null})
    var histoFiles = Array<File?>(sampleNames.size,{null})
    val tmp = File(Desirable.tmpDir,"KMCtmp")

    init {sampleNames = compactSampleNames(sampleNames)}

    private fun compactSampleNames(sampleNames: Array<String>): Array<String> {
        val minLen = sampleNames.map { it.length }.min()!!
        val start = (0 until minLen).find{ i ->
            val column = sampleNames.map { it[i] }
            column.any { it!=column[0] }
        }!!
        val end = (1..minLen).find{ i ->
            val column = sampleNames.map { it[it.length-i] }
            column.any { it!=column[0] }
        }!!
        return sampleNames.map { it.substring(start).dropLast(end-1) }.toTypedArray()
    }

    private fun runKmerCount(): Array<File?> =
        readFiles.mapIndexed { i, sample ->
            val sampleFile = File(Desirable.tmpDir,sampleNames[i])
            sampleFile.writeText(sample.joinToString("\n") { file->file.absolutePath })
            val cmd = "${Desirable.KMCPath} -k${Desirable.K} -cs1000000000 @$sampleFile $sampleFile $tmp"
            Desirable.logRuntime("KMC","cmd: $cmd\nCount Kmer: ${sampleNames[i]}",cmd)
            sampleFile
        }.toTypedArray()

    private fun runKmerHisto(): Array<File?> =
        countFiles.mapIndexed { i, sample->
            val histoFile = File(Desirable.tmpDir,sampleNames[i]+".histo")
            val cmd = "${Desirable.KMCToolsPath} transform $sample histogram $histoFile -cx1000000"
            Desirable.logRuntime("KMC","cmd: $cmd\nKmer Histogram: ${sampleNames[i]}",cmd)
            histoFile
        }.toTypedArray()

    fun getSpectrumTroughPeak(histoFile: File): Pair<Int, Int> { // get first local minimum in kmer spectrum
        val spectrum = histoFile.readLines().map {
            val tokens = it.split('\t');
            Pair(tokens.first().toLong(),tokens.last().toLong())
        }.toTypedArray()
        val total = spectrum.map { it.first * it.second }.sum()
        var accumulate = 0L
        var trough = spectrum.indices.first {
            if (spectrum[it].second <= spectrum[it + 1].second) return@first true
            accumulate += spectrum[it].first * spectrum[it].second
            (accumulate > total * 0.5)
        }
        if (accumulate > total * 0.5) trough = 0
        val peak = spectrum.filter { it.first>=trough }.maxBy { it.first * it.second }!!.first
        println("[KmerOperator]\n In ${histoFile.name}:\nTrough (Error threshold): $trough\nPeak (Coverage): $peak")
        return Pair(trough.toInt(), peak.toInt())
    }

    private fun getExclusiveKmerInTarget(countFiles: Array<File?>, thoughPeaks: Array<Pair<Int, Int>>): File {
        val outputFile = File(Desirable.tmpDir,"ExclusiveKmer")
        val complexScript = File(Desirable.tmpDir,"ComplexScript.kmc")
        val inputList = countFiles.mapIndexed { i, file -> "set${i+1} = $file -ci${thoughPeaks[i].first}" }.joinToString(separator = "\n")
        val outputOps = "$outputFile = ${countFiles.indices.joinToString(separator = " - ") { "set${it+1}" }}"
        complexScript.writeText("INPUT:\n$inputList\nOUTPUT:\n$outputOps\nOUTPUT_PARAMS:\n")
        var cmd = "${Desirable.KMCToolsPath} complex $complexScript"
        Desirable.logRuntime("KMC","cmd: $cmd\nScreen Exclusive Kmer in ${sampleNames[0]}",cmd)
        cmd = "${Desirable.KMCToolsPath} transform $outputFile dump $outputFile"
        Desirable.logRuntime("KMC","cmd: $cmd\nDump Kmer in ${sampleNames[0]}",cmd)
        return outputFile
    }

    private fun getExclusiveReadInTarget(countFiles: File, exKmerFile: File): File {
        val outputFile = File(Desirable.tmpDir,"ExclusiveRead.fastq")
        val cmd = "${Desirable.KMCToolsPath} filter $exKmerFile @$countFiles -ci0.8 $outputFile"
        Desirable.logRuntime("KMC","cmd: $cmd\nScreen Exclusive Read in ${sampleNames[0]}",cmd)
        return outputFile
    }

    private fun getExtendedKmerInReads(exReadFile: File, countFile: File?, ci:Int): File {
        val rawKmer = File(Desirable.tmpDir,"rawKmer")
        var cmd = "${Desirable.KMCPath} -k${Desirable.K} -ci1 -cs1000000000 $exReadFile $rawKmer $tmp"
        Desirable.logRuntime("KMC","cmd: $cmd\nCount raw extended Kmer in ${sampleNames[0]}",cmd)
        val outputFile = File(Desirable.tmpDir,"ExtendedKmer")
        cmd = "${Desirable.KMCToolsPath} simple $rawKmer $countFile intersect $outputFile -ocmax -ci$ci"
        Desirable.logRuntime("KMC","cmd: $cmd\nScreen solid extended Kmer in ${sampleNames[0]}",cmd)
        cmd = "${Desirable.KMCToolsPath} transform $outputFile dump $outputFile"
        Desirable.logRuntime("KMC","cmd: $cmd\nDump Kmer in ${sampleNames[0]}",cmd)
        return outputFile
    }

    fun run(): List<File> {
        if (!tmp.exists()) tmp.mkdir()
        countFiles = runKmerCount()
        histoFiles = runKmerHisto()
        val thoughPeaks = histoFiles.map { getSpectrumTroughPeak(it!!) }.toTypedArray()
        val exKmerFile = getExclusiveKmerInTarget(countFiles,thoughPeaks)
        val exReadFile = getExclusiveReadInTarget(countFiles[0]!!,exKmerFile)
        val extendedKmer = getExtendedKmerInReads(exReadFile, countFiles[0]!!,thoughPeaks[0].first)
        return listOf(extendedKmer,exReadFile)
    }

}
