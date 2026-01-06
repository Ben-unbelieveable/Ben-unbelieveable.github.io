version 1.0
## language:   WDL
## File    :   main.wdl
## Time    :   2024/05/09 13:40:01
## Author  :   Liu.Bo
## Version :   v1.0
## Contact :   liubo4@genomics.cn/614347533@qq.com
## WebSite :   http://www.ben-air.cn/
##
## Cromwell version support :
##            - Successfully tested on v36
##            - Does not work on versions < v23 due to output syntax
##
## WORKFLOW DEFINITIONS
# description

workflow Main_workflow {
 call docker_test {
  input:
 }
 output{
   File OUTA = docker_test.output_result
 }
}


task docker_test {
  command <<<
    ls / > list
  >>>
  runtime {
  docker: "centos"
  }
  output {
    File output_result = "list"
  }
}
