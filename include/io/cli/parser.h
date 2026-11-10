//
// Created by Mert on 30.10.26.
//

#ifndef E_ROUTING_PARSER_H
#define E_ROUTING_PARSER_H


#include "io/cli/command.h"
#include "core/errors.h"

Result<Command> parse(int argc, char** argv);


#endif // E_ROUTING_PARSER_H
