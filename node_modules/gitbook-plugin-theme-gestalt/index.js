module.exports =
{
    hooks:
    {
        config: function(config)
        {
            config.styles = config.styles || config.pluginsConfig['theme-default'].styles;
            console.log( config.styles );
            return config;
        }
    }
};


