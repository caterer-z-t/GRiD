version 1.0

workflow HelloInput {
    input {
        String name
    }


  call WriteGreeting { input: name = name }
}

task WriteGreeting {
  input {
    String name
  }
  command {
    echo "hello ${name}!"
  }
  output {
    File response = stdout()
  }
}
