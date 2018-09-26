package edu.nus

import com.github.ajalt.clikt.core.*
import com.github.ajalt.clikt.parameters.arguments.argument
import com.github.ajalt.clikt.parameters.arguments.multiple
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.*
import edu.gatech.kanalyze.module.KAnalyzeModule
import java.io.File

class Desirable : CliktCommand() {
    val threads by option(help = "Threads to use, default by the CPU core number").int().default(Runtime.getRuntime().availableProcessors())
    val K by option("-k", help = "Kmer size").int().default(31).validate { if (it>32) fail("only support K <= 32") else Util.K = it }
    val output by option("-o", help = "Output file path").file().default(File("DesirableOut.fasta"))
    val tmpDir: File by option("-d", help = "Temporary folder for intermediate results").file(exists = false).default(File("DesirableTmp"))
    val target: List<File> by option("-t", help = "Mark next file is target sample").file(exists = true).multiple().validate { it.isNotEmpty() }
    val control: List<File> by argument().file(exists = true).multiple()
    var targetCount: File? = null
    var controlCount: File? = null

    fun countKmerInFile(files: List<File>):File {
        val filename = files.first().name
        val countFile = File(tmpDir, "$filename.kc")
        val format = (if (filename.contains("fastq") || filename.contains("fq")) "fastq" else "fasta") + if (filename.endsWith(".gz")) "gz" else ""
        val arguments = "count -rcanonical " +
                "-t $threads " +
                "-k $K " +
                "-c kmercount:2 " +
                "-o ${countFile.path} " +
                "-f $format ${files.map { it.path } .joinToString(separator=" ")}"
        println("[KmerCount] start\n Command: $arguments")
        val runtime = kotlin.system.measureTimeMillis { KAnalyzeModule.main(arguments.split(" ").toTypedArray()) }
        println("[KmerCount] finished in ${runtime/1000} seconds")
        return countFile
    }

    override fun run() {
        tmpDir.mkdirs()
        targetCount = countKmerInFile(target)
        controlCount = countKmerInFile(control)
        val ka = KmerAnalyzer(this)
        val targetTroughPeak = ka.getSpectrumTroughPeak(targetCount!!)
        val controlTroughPeak = ka.getSpectrumTroughPeak(controlCount!!)
        val exclusiveKmers = ka.filterTargetByControl(target, targetTroughPeak, controlTroughPeak)
        val assembler = Assembler(this)
        assembler.constructGraphFromKmerSet(exclusiveKmers)
        assembler.getContigByTraverse()
        assembler.writeContigsToFasta(200)
//        tmpDir.deleteRecursively()
    }
}

fun main(args: Array<String>) = Desirable().main(args)