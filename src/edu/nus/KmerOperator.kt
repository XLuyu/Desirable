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
        var end = (1..minLen).find{ i ->
            val column = sampleNames.map { it[it.length-i] }
            column.any { it!=column[0] }
        }!!
        if (start+end-1>=minLen) end = 1
        return sampleNames.map { it.substring(start).dropLast(end-1) }.toTypedArray()
    }

    private fun runKmerCount(): Array<File?> =
        readFiles.mapIndexed { i, sample ->
            val sampleFile = File(Desirable.tmpDir,sampleNames[i])
            sampleFile.writeText(sample.joinToString("\n") { file->file.absolutePath })
            val cmd = "${Desirable.KMCPath} -k${Desirable.K} -cs1000000000 @$sampleFile $sampleFile $tmp"
            Desirable.logRuntime("KMC/Count","<${sampleNames[i]}> \"$cmd\"\n",cmd)
            sampleFile
        }.toTypedArray()

    private fun runKmerHisto(): Array<File?> =
        countFiles.mapIndexed { i, sample->
            val histoFile = File(Desirable.tmpDir,sampleNames[i]+".histo")
            val cmd = "${Desirable.KMCToolsPath} transform $sample histogram $histoFile -cx1000000"
            Desirable.logRuntime("KMC/Histogram","<${sampleNames[i]}> \"$cmd\"\n",cmd)
            histoFile
        }.toTypedArray()

    fun getSpectrumTroughPeak(histoFile: File): Pair<Int, Int> { // get first local minimum in kmer spectrum
        val spectrum = histoFile.readLines().map {
            line -> line.split('\t').let { Pair(it.first().toLong(),it.last().toLong()) }
        }.filter { it.second!=0L } .toTypedArray()
        spectrum.sortBy { it.first }
        val total = spectrum.map { it.first * it.second }.sum()
        var accumulate = 0L
        val N70 = spectrum.first {
            accumulate += it.first * it.second
            (accumulate > total * 0.3)
        }
        var trough = spectrum.filter { it.first<=N70.first }.minBy { it.second }!!.first
        if (trough == N70.first) trough = 0
        val peak = spectrum.filter { it.first>trough }.maxBy { it.first * it.second }!!.first
        println("[KmerOperator] <${histoFile.name}>\nKmer N70: $N70\nTrough (Error threshold): $trough\nPeak (Coverage): $peak")
        return Pair(trough.toInt(), peak.toInt())
    }

    private fun getExclusiveKmerInTarget(countFiles: Array<File?>, thoughPeaks: Array<Pair<Int, Int>>): File {
        val outputFile = File(Desirable.tmpDir,"ExclusiveKmer")
        val complexScript = File(Desirable.tmpDir,"ComplexScript.kmc")
        val inputList = countFiles.mapIndexed { i, file -> "set${i+1} = $file -ci${thoughPeaks[i].first}" }.joinToString(separator = "\n")
        val outputOps = "$outputFile = ${countFiles.indices.joinToString(separator = " - ") { "set${it+1}" }}"
        complexScript.writeText("INPUT:\n$inputList\nOUTPUT:\n$outputOps\n")
        var cmd = "${Desirable.KMCToolsPath} complex $complexScript" // exclude positive solid kmer by negative solid kmer
        Desirable.logRuntime("KMC/Complex","<${sampleNames[0]}> Exclusive Kmer \"$cmd\"\n",cmd)
        cmd = "${Desirable.KMCToolsPath} transform $outputFile dump $outputFile"
        Desirable.logRuntime("KMC/Dump","<${sampleNames[0]}> \"$cmd\"\n",cmd)
        if (outputFile.length()<Desirable.K) {
            println("[KMC] No exclusive K-mer. Please try a larger K.")
            System.exit(1)
        }
        return outputFile
    }

    private fun getExclusiveReadInTarget(countFiles: File, exKmerFile: File): File {
        val outputFile = File(Desirable.tmpDir,"ExclusiveRead.fastq")
        val cmd = "${Desirable.KMCToolsPath} filter $exKmerFile @$countFiles -ci0.8 $outputFile" // fetch reads with >80% its kmers are exclusive, @countFiles means a list
        Desirable.logRuntime("KMC/FilterRead","<${sampleNames[0]}> Exclusive Read \"$cmd\"\n",cmd)
        if (outputFile.length()<Desirable.K) {
            println("[KMC] No exclusive sequencing as exclusive k-mers are dispersed. Please try a larger K.")
            System.exit(1)
        }
        return outputFile
    }

    private fun getExtendedKmerInReads(exReadFile: File, countFile: File?, ci:Int): File {
        val rawKmer = File(Desirable.tmpDir,"rawKmer")
        var cmd = "${Desirable.KMCPath} -k${Desirable.K} -ci1 -cs1000000000 $exReadFile $rawKmer $tmp" // KMC the subset reads
        Desirable.logRuntime("KMC/Count","<${sampleNames[0]}> \"$cmd\"\n",cmd)
        val outputFile = File(Desirable.tmpDir,"ExtendedKmer")
        cmd = "${Desirable.KMCToolsPath} simple $rawKmer $countFile intersect $outputFile -ocmax -ci$ci" // filter out solid kmer
        Desirable.logRuntime("KMC/Intersect","<${sampleNames[0]}> Solid Extended Kmer \"$cmd\"\n",cmd)
        cmd = "${Desirable.KMCToolsPath} transform $outputFile dump $outputFile"
        Desirable.logRuntime("KMC/Dump","<${sampleNames[0]}> \"$cmd\"\n",cmd)
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
