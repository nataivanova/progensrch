#!/usr/bin/env python3

import os
import uuid
import logging
from progenlib import *
from telegram.ext import Updater, MessageHandler, CommandHandler, Filters

logger = logging.getLogger('progen tool - telegram')
req_id = 'startup'
logging.basicConfig( level=logging.INFO
                    , format = 'id=' +  req_id + ', t=%(asctime)s, level=%(levelname)s, message=%(message)s' )

BOT_TOKEN = os.environ.get('TELEGRAM_TOKEN')
db_location = os.environ['DB_LOCATION']

db = ProgenitorDatabase(db_location)

def downloader(update, context):

    req_id = str(uuid.uuid4())
    logging.basicConfig( level=logging.INFO
                    , format = 'id=' +  req_id + ', t=%(asctime)s, level=%(levelname)s, message=%(message)s' )

    logger.info('received a file, addigning request id=' + req_id )

    reqfile = '/tmp/' + req_id
    with open(reqfile, 'wb') as _:
        context.bot.get_file(update.message.document).download( out = _ )

    context.bot.send_message(
        chat_id = update.message.chat_id
        , text    = "your search id is " + req_id)
    try:
        query = ProgenitorQuery(reqfile)
        results = ProgenitorSearch(query, db)
        outpath = '/tmp/' + req_id + '.result'
        try:
            logger.debug('saving search results to a file for sending')
            outfile = open(outpath, 'w')
            outfile.write( str(results) )
            outfile.close()
            outfile = open(outpath, 'rb')
            context.bot.send_document(
                chat_id  = update.message.chat_id
                , document = outfile
                , disable_content_type_detection = True)
        except Exception as e:
            logger.error ('failed to create output file')
            logger.error (str(e))

    except Exception as e:
        logger.error ('could not parse query or execute search')
        logger.error (str(e))
        try:
            errpath = '/tmp/' + req_id + '.err'
            errfile = open(errpath, 'wb')
            errfile.write(bytes(str(e),'UTF-8'))
            errfile.close()
            errfile = open(errpath, 'rb')
            context.bot.send_document(
                chat_id = update.message.chat_id
                , document=errfile
                , disable_content_type_detection=True)
        except Exception as e:
            logger.critical ('Failed to write error to file')
            logger.critical(str(e))

def help(update, context):
    logger.info("help requested")
    x = "help will be given to those who ask for it"
    context.bot.send_message( chat_id = update.message.chat_id
                              , text = x)
    try:
        with open ('./sample_query.txt', 'rb') as sample_file:
            context.bot.send_document ( chat_id = update.message.chat_id
                                        , document = sample_file
                                        , disable_content_type_detection = True)
    except Exception as e:
        logger.error ('cannot send sample query file')
        logger.error (str(e))

updater = Updater(BOT_TOKEN, use_context=True)

updater.dispatcher.add_handler(MessageHandler(Filters.document, downloader))
updater.dispatcher.add_handler(CommandHandler("help", help))

updater.start_polling()
updater.idle()
