plugins {
    // using Java plugin instead of Java application
    java
    // since we take user input via stdin
    // recommend to build and directly run the jar rather than using a Gradle application wrapper
    // `gradle jar`
    // `java -jar ./app/build/libs/app.jar`
    jacoco
}
repositories {
    // Use Maven Central for resolving dependencies.
    mavenCentral()
}

dependencies {
    // Use JUnit Jupiter for testing.
    testImplementation(libs.junit.jupiter)

    testRuntimeOnly("org.junit.platform:junit-platform-launcher")

    // This dependency is used by the application.
    implementation(libs.guava)

}

// Apply a specific Java toolchain to ease working on different environments.
java {
    toolchain {
        languageVersion = JavaLanguageVersion.of(21)
    }
}

tasks.jar {
  //https://docs.gradle.org/9.0.0/dsl/org.gradle.jvm.tasks.Jar.html#org.gradle.jvm.tasks.Jar:manifest
  //https://docs.gradle.org/9.0.0/userguide/building_java_projects.html#sec:jar_manifest
  manifest {
    attributes("Main-Class" to "org.wbdbga.DeBruijnGenomeAssembler")
  }
}

tasks.named<Test>("test") {
    // Use JUnit Platform for unit tests.
    useJUnitPlatform()
}

tasks.test {
    finalizedBy(tasks.jacocoTestReport) // report is always generated after tests run
}
tasks.jacocoTestReport {
    dependsOn(tasks.test) // tests are required to run before generating the report
}