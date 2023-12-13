# Writing Modules

\brief Instructions for writing custom ampsci modules

[[Home](/README.md)]

- The modules system allows the easy calculation of any atomic properties after the wavefunction has been calculated.
- Any number of _modules_ can be run by adding a `Module::moduleName{}' block to the input file.
- Get a list of available modules: `./ampsci -m`
- See [doc/modules.md](/doc/modules.md) for details of currently available modules
- The code is designed so that you can easily create your own modules.

## Creating your own module

An example module is provided to help you write your own module;

- [src/Modules/exampleModule.hpp](/src/Modules/exampleModule.hpp)
- [src/Modules/exampleModule.cpp](/src/Modules/exampleModule.cpp)

You should duplicate this module (both the .cpp and .hpp files) and give it a new name. That will be much easier than starting from scratch.

- Modules are functions that have the following function signature:

  ```cpp
  namespace Module{
  void exampleModule(const IO::InputBlock &input, const Wavefunction &wf);
  }
  ```

- `input` is an `IO::InputBlock` that holds any input options.

- `wf` is the `Wavefunction` object that was calculated by ampsci

- They are typically placed in the `Module` namespace (but don't need to be).

- Typically, modules live in the `src/Modules` directory, but they can live anywhere

- Then, you can do anything you like inside this function

## Including your Module into ampsci

- In order for ampsci to know about your module, you must update the file `src/Modules/module_list.hpp`
  - Add the corresponding header file to the #includes list

- Add a `std::pair` the the `std::vector` _module_list_ in the form

    ```cpp
    {"moduleName", &moduleName}
    ```

- The first of the pair is a string, which will be the name of the module. This is how you will refer to the module in the input file.

- The second of the pair is the pointer to the function name.

- You then have to recompile ampsci (modules are compiled into ampsci for now)

- That's it - you're now ready to run your module by adding the `Module::moduleName{}` block to the input file.

## Highly recommended (but optional)

- It's highly recommended that you add a 'check()' statement for any input options that you use in your module (see example below)
- This has two benefits:
  - Firstly, it checks for possible spelling mistakes in user inputs
  - (If an option) is spelled incorrectly, it will otherwise be ignored. This eaves the user thinking they set an option when they haven't
  - Secondly, it allows you to provide a short description of each option, which will be printed to the screen when the user requests 'help' for a given Module

```cpp
  // Check the input option for spelling mistakes + provide description
  input.check({{"option1", "Short description of option1 [default]"},
               {"option2", "Short description of option2 [default]"}});
```

- Finally, it's strongly recommended to immediately exit the Module after the check() if 'help' was requested (see example below). This just limits unwanted noise/screen output
- i.e., we generally don't want to actually run the module if we were just requesting help for it

```cpp
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }
```

- Minimal example:

```cpp
  void exampleModule(const IO::InputBlock &input, const Wavefunction &wf){

    input.check({{"option1", "Short description of option1 [default1]"},
               {"option2", "Short description of option2 [default2]"}});

    // If we are just requesting 'help', don't run module:
    if (input.has_option("help")) {
      return;
    }

    // read in the input options:
    auto option1 = input.get("option1", default1);
    auto option2 = input.get("option2", default2);
    // The variable 'option1' will be set according to the user input
    // If no user input for option1 was given, it will be set to default1
    // See documentation for IO::InputBlock for more detail
  }
```
